#
#  farPairs.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 7/13/10.
#

try:
    import psyco
    psyco.full()
except:
    pass

import sys, time
import optparse
from commoncode import readDataset

print "%prog: version 1.3"


def main(argv=None):
    if not argv:
        argv = sys.argv

    usage = "usage: python %prog rdsfile outfile bedfile [--verbose] [--cache numPages] [--minDist bp] [--maxDist bp] [--minCount count] [--label string]"

    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--sameChromOnly", action="store_true", dest="sameChromOnly")
    parser.add_option("--cache", type="int", dest="cachePages")
    parser.add_option("--verbose", action="store_true", dest="doVerbose")
    parser.add_option("--minDist", type="int", dest="minDist")
    parser.add_option("--maxDist", type="int", dest="maxDist")
    parser.add_option("--minCount", type="int", dest="minCount")
    parser.add_option("--label", dest="label")
    parser.set_defaults(sameChromOnly=False, doVerbose=False, cachePages=None,
                        minDist=1000, maxDist=500000, minCount=2, label=None)
    (options, args) = parser.parse_args(argv[1:])

    if len(args) < 3:
        print usage
        print "\tIs both slow and takes up large amount of RAM"
        sys.exit(1)

    rdsfile = args[0]
    outfilename = args[1]
    outbedname = args[2]

    farPairs(rdsfile, outfilename, outbedname, options.sameChromOnly, options.doVerbose,
             options.cachePages, options.minDist, options.maxDist, options.minCount,
             options.label)


def farPairs(rdsfile, outfilename, outbedname, sameChromOnly=False, doVerbose=False,
             cachePages=None, minDist=1000, maxDist=500000, minCount=2, label=None):

    doCache = False
    if cachePages is not None:
        doCache = True
    else:
        cachePages = 0

    if label is None:
        label = rdsfile

    RDS = readDataset(rdsfile, verbose=True, cache=doCache)
    rdsChromList = RDS.getChromosomes()

    if doVerbose:
        print time.ctime()

    total = 0
    outfile = open(outfilename, "w")
    outbed = open(outbedname, "w")
    outbed.write('track name="%s distal pairs" color=0,255,0\n' % label)

    readlen = RDS.getReadSize()
    flagDict = {}
    for chromosome in rdsChromList:
        if doNotProcessChromosome(chromosome):
            continue

        print chromosome
        uniqDict = RDS.getReadsDict(fullChrom=True, chrom=chromosome, noSense=True, withFlag=True, withPairID=True, doUniqs=True, readIDDict=True)
        if doVerbose:
            print len(uniqDict), time.ctime()

        for readID in uniqDict:
            readList = uniqDict[readID]
            if len(readList) == 2:
                total += 1
                (start1, flag1, pair1) = readList[0]
                (start2, flag2, pair2) = readList[1]

                if flag1 != flag2:
                    dist = abs(start1 - start2)
                    startList = [start1, start2]
                    stopList = [start1 + readlen, start2 + readlen]
                    startList.sort()
                    stopList.sort()
                    if flag1 != "" and flag2 != "" and minDist < dist < maxDist:
                        outputLine = splitReadWrite(chromosome, 2, startList, stopList, "+", readID, "0,255,0", "0,255,0")
                        outbed.write(outputLine)
                        if doVerbose:
                            print flag1, flag2, dist

                        try:
                            flagDict[flag1].append((flag2, start1, start2))
                        except KeyError:
                            flagDict[flag1] = [(flag2, start1, start2)]

                        try:
                            flagDict[flag2].append((flag1, start1, start2))
                        except KeyError:
                            flagDict[flag2] = [(flag2, start1, start2)]

    print "%d connected regions" % len(flagDict)

    for region in flagDict:
        flagDict[region].sort()
        regionConnections = {}
        for (region2, start1, start2) in flagDict[region]:
            try:
                regionConnections[region2] += 1
            except KeyError:
                regionConnections[region2] = 1

        for region2 in regionConnections:
            if regionConnections[region2] >= minCount:
                outfile.write("%s\t%s\t%d\n" % (region, region2, regionConnections[region2]))
                if doVerbose:
                    print "%s\t%s\t%d" % (region, region2, regionConnections[region2])

    outfile.close()
    outbed.close()
    if doVerbose:
        print "finished: ", time.ctime()


def doNotProcessChromosome(chrom):
    return chrom == "chrM"


def splitReadWrite(chrom, numPieces, startList, stopList, rsense, readName, plusSense, minusSense):
    readSizes = "%d" % (stopList[0] - startList[0])
    readCoords = "0"
    leftStart = startList[0] - 1
    rightStop = stopList[-1]
    for index in range(1, numPieces):
        readSizes += ",%d" % (stopList[index] - startList[index] + 1)
        readCoords += ",%d" % (startList[index] - startList[0])

    if rsense == "+":
        senseCode = plusSense
    else:
        senseCode = minusSense

    outline = "%s\t%d\t%d\t%s\t1000\t%s\t0\t0\t%s\t%d\t%s\t%s\n" % (chrom, leftStart, rightStop, readName, rsense, senseCode, numPieces, readSizes, readCoords)
    return outline


if __name__ == "__main__":
    main(sys.argv)