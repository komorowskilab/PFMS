#
#  distalPairs.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 10/14/08.
#


try:
    import psyco
    psyco.full()
except:
    pass

from commoncode import *
import sys, time

print '%s: version 3.3' % sys.argv[0]

if len(sys.argv) < 2:
    print 'usage: python %s minDist rdsfile outfile [-sameChrom] [-splices] [-maxDist bp] [-verbose] [-cache cachepages]' % sys.argv[0]
    print '\tlooks at all chromosomes simultaneously: is both slow and takes up large amount of RAM'
    sys.exit(1)

minDist = int(sys.argv[1])
rdsfile = sys.argv[2]
outfilename = sys.argv[3]

sameChromOnly = False
if '-sameChrom' in sys.argv:
    sameChromOnly = True

doSplices = False
if '-splices' in sys.argv:
    doSplices = True
    
doVerbose = False
if '-verbose' in sys.argv:
    doVerbose = True

maxDist = 1000000000
if '-maxDist' in sys.argv:
    maxDist = int(sys.argv[sys.argv.index('-maxDist') + 1])

doCache = False
cachePages = -1
if '-cache' in sys.argv:
    doCache = True
    try:
        cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])
    except: 
        pass

RDS = readDataset(rdsfile, verbose = True, cache=doCache)
if not RDS.hasIndex():
    print "Will not attempt to run on unIndexed dataset - please index with rdsmetadata.py and rerun"
    sys.exit(1)


if cachePages > RDS.getDefaultCacheSize():
    RDS.setDBcache(cachePages)

print time.ctime()

if doSplices:
    print "getting splices"
    splicesDict = RDS.getSplicesDict(withChrom=True, withPairID=True, readIDDict=True, splitRead=True)
    print "got splices"

print "getting uniq reads"    
uniqDict = RDS.getReadsDict(withChrom=True, withPairID=True, doUniqs=True, readIDDict=True)
print "got uniqs"

if doSplices:
    for readID in splicesDict:
        theRead = splicesDict[readID]
        read0 = theRead[0]
        del read0[1]
        try:
            uniqDict[readID].append(read0)
        except:
            if len(theRead) == 4:
                read2 = theRead[2]
                del read2[1]
                uniqDict[readID] = [read0,read2]
if doVerbose:
    print len(uniqDict), time.ctime()

outfile = open(outfilename,'w')

diffChrom = 0
distal = 0
total = 0
for readID in uniqDict:
    readList = uniqDict[readID]
    if len(readList) == 2:
        total += 1
        (start1, sense1, chrom1, pair1) = readList[0]
        (start2, sense2, chrom2, pair2) = readList[1]
        
        if chrom1 != chrom2:
            diffChrom += 1
            if sameChromOnly:
                continue
            else:
                outline = '%s\t%s\t%d\t%s\t%s\t%d\t%s' % (readID, chrom1, start1, sense1, chrom2, start2, sense2)
                outfile.write(outline + '\n')
                if doVerbose:
                    print diffChrom, outline
        else:
            dist = abs(start1 - start2)
        
            if minDist < dist < maxDist:
                distal += 1
                outline = '%s\t%s\t%d\t%s\t%d\t%s\t%d' % (readID, chrom1, start1, sense1, start2, sense2, dist)
                outfile.write(outline + '\n')
                if doVerbose:
                    print distal, outline

outfile.write('#distal: %d\tdiffChrom: %d\tpossible: %d\n' % (distal, diffChrom, total))
total = float(total)
if total < 1:
    total = 1.
outfile.write('#distal %2.2f pct\tdiffChrom %2.2f pct\n' % ((100. * distal/total), (100. * diffChrom/total)))
outfile.close()
print "distal: %d\tdiffChrom: %d\tpossible: %d" % (distal, diffChrom, int(total))
print "distal: %2.2f pct\tdiffChrom: %2.2f pct\n" % ((100. * distal/total), (100. * diffChrom/total))
print time.ctime()


