"""
MakeBamFromRds

Converts ERANGE RDS zero based file to Bam zero based format.

Usage: python MakeBamFromRDS.py rdsFile bamFile [options]

"""

try:
    import psyco
    psyco.full()
except:
    pass

import sys, re, string, optparse, random
import pysam
from commoncode import readDataset


def main(argv=None):
    if not argv:
        argv = sys.argv

    verstring = "%prog: version 3.1"
    print verstring

    doPairs = False
    
    usage = "usage: python %prog rdsFile bamFile [options]"

    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--nouniq", action="store_false", dest="withUniqs")
    parser.add_option("--nomulti", action="store_false", dest="withMulti")
    parser.add_option("--splices", action="store_true", dest="doSplices")
    parser.add_option("--flag", dest="withFlag")
    parser.add_option("--flaglike", action="store_true", dest="useFlagLike")
    parser.add_option("--pairs", action="store_true", dest="doPairs")
    parser.add_option("--cache", type="int", dest="cachePages")
    parser.add_option("--enforceChr", action="store_true", dest="enforceChr")
    parser.add_option("--chrom", action="append", dest="chromList")
    parser.add_option("--fasta", dest="fastaFileName")
    parser.set_defaults(withUniqs=True, withMulti=True, doSplices=False,
                        doPairs=False, withFlag="", useFlagLike=False, enforceChr=False,
                        doCache=False, cachePages=100000, fastaFileName="",
                        chromList=[])

    (options, args) = parser.parse_args(argv[1:])

    if len(args) < 2:
        print usage
        sys.exit(1)

    rdsfile = args[0]
    outfilename = args[1]

    allChrom = True
    if options.chromList:
        allChrom = False

    makeBamFromRds(rdsfile, outfilename, options.withUniqs, options.withMulti,
                     options.doSplices, doPairs, options.withFlag, options.useFlagLike,
                     options.enforceChr, allChrom, options.doCache, options.cachePages,
                     options.chromList, options.fastaFileName)


def makeBamFromRds(rdsfile, outfilename, withUniqs=True, withMulti=True,
                     doSplices=False, doPairs=False, withFlag="",
                     useFlagLike=False, enforceChr=False, allChrom=True,
                     doCache=False, cachePages=100000, chromList=[], fastaFileName=""):

    if not withUniqs and not withMulti and not doSplices:
        print "must be outputting at least one of uniqs, multi, or -splices - exiting"
        sys.exit(1)

    print "\nsample:"
    RDS = readDataset(rdsfile, verbose = True, cache=doCache)

    if cachePages > RDS.getDefaultCacheSize():
        RDS.setDBcache(cachePages)

    readlength = RDS.getReadSize()

    if allChrom:
        if withUniqs:
            chromList = RDS.getChromosomes()
        elif withMulti:
            chromList = RDS.getChromosomes(table="multi")
        else:
            chromList = RDS.getChromosomes(table="splices")

        chromList.sort()

    fastaSequenceDict = {}
    if fastaFileName:
        fastafile = open(fastaFileName)
        fastaSequenceDict = getFastaSequenceDictionary(fastaFileName)
        fastafile.close()

    referenceSequenceList = []
    chromRemoveList = []
    for chromosome in chromList:
        if doNotOutputChromosome(chromosome, enforceChr):
            chromRemoveList.append(chromosome)
        else:
            chromosomeLength = RDS.getMaxCoordinate(chromosome, doUniqs=withUniqs, doMulti=withMulti, doSplices=doSplices)
            referenceDataDict = {"LN": int(chromosomeLength), "SN": str(chromosome)}
            referenceSequenceList.append(referenceDataDict)

    for chrom in chromRemoveList:
        chromList.remove(chrom)

    header = {"HD": {"VN": "1.0"}}
    if referenceSequenceList:
        header["SQ"] = referenceSequenceList

    outfile = pysam.Samfile(outfilename, "wb", header=header)

    totalWrites = 0
    noncanonicalSplices = 0
    for chrom in chromList:
        index = 0
        print "chromosome %s" % (chrom)
        if withUniqs or withMulti:
            hitDict = RDS.getReadsDict(fullChrom=True, chrom=chrom, flag=withFlag, withWeight=True, withID=True,
                                       withPairID=doPairs, doUniqs=withUniqs, doMulti=withMulti, readIDDict=False,
                                       flagLike=useFlagLike, entryDict=True)

            for read in hitDict[chrom]:
                index += writeBAMEntry(outfile, chrom, read, readlength)

        if doSplices:
            numSpliceReadsWritten, noncanonical = processSpliceReads(RDS, outfile, chrom, withFlag, useFlagLike, readlength, fastaSequenceDict)
            index += numSpliceReadsWritten
            noncanonicalSplices += noncanonical

        print index
        totalWrites += index

    outfile.close()
    print "%d total reads written" % totalWrites
    print "%d non-canonical splices" % noncanonicalSplices


def processSpliceReads(RDS, outfile, chrom, withFlag, useFlagLike, readlength, fastaSequenceDict={}):
    index = 0
    noncanonicalSplices = 0
    spliceDict = RDS.getSplicesDict(fullChrom=True, chrom=chrom, flag=withFlag, withID=True, flagLike=useFlagLike, entryDict=True, withWeight=True)
    if chrom not in spliceDict:
        pass
    else:
        for read in spliceDict[chrom]:
            if fastaSequenceDict.has_key(chrom):
                read["sense"], noncanonical = fixSpliceSense(fastaSequenceDict[chrom], chrom, read["startR"], read["stopL"], read["sense"])
                noncanonicalSplices += noncanonical

            index += writeBAMEntry(outfile, chrom, read, readlength)

    return index, noncanonicalSplices


def writeBAMEntry(outfile, chrom, outputDict, readlength):
    index = 0
    tagList = []
    alignedRead = pysam.AlignedRead()
    queryName = string.split(outputDict["readID"], "/")[0]
    alignedRead.qname = queryName
    if outputDict["sense"] == "-":
        alignedRead.is_reverse = True

    alignedRead.rname = outfile.references.index(chrom)

    if outputDict.has_key("startL"):
        startL = outputDict["startL"]
        stopL = outputDict["stopL"]
        startR = outputDict["startR"]
        stopR = outputDict["stopR"]
        alignedRead.pos = startL
        alignedRead.cigar = [(0,stopL - startL + 1), (3, startR - stopL - 1), (0, stopR - startR + 1)]
        tagList.append(("XS", outputDict["sense"]))
    else:
        alignedRead.pos = outputDict["start"]
        alignedRead.cigar = [(0, readlength)]

    if outputDict.has_key("pairID"):
        pairID = outputDict["pairID"]
        if pairID == "1":
            alignedRead.is_read1 = True
            alignedRead.is_proper_pair = True
        elif pairID == "2":
            alignedRead.is_read2 = True
            alignedRead.is_proper_pair = True
        else:
            pass

    if outputDict.has_key("mismatch"):
        mismatchTag = getMismatches(outputDict["mismatch"])
        if mismatchTag:
            tagList.append(("MD", mismatchTag))
    
    if tagList:
        alignedRead.tags = tagList

    multiplicity = 1.0 / outputDict.get("weight", 1.0)
    while multiplicity > 0:
        outfile.write(alignedRead)
        multiplicity -= 1.0
        index += 1

    return index


def getMismatches(mismatchString):
    mismatch = ""
    positions = re.findall("\d+", mismatchString)
    nucleotides = re.findall("([ACGTN])\d+", mismatchString)
    for index in range(0, len(positions)):
        mismatch = "%s%s%s" % (mismatch, positions[index], nucleotides[index])

    return mismatch


def doNotOutputChromosome(chrom, enforceChr):
    result = False

    if chrom == "chrM":
        result = True

    if enforceChr and ("chr" not in chrom):
        result = True

    return result


def getFastaSequenceDictionary(fastaFileName):
    fastaSeqDict = {}
    fchrom = ""
    fseq = ""

    fastafile = open(fastaFileName)
    for line in fastafile:
        if line[0] == ">":
            if fchrom != "":
                fastaSeqDict[fchrom] = fseq

            fseq = ""
            fchrom = line[1:-1]
        else:
            fseq += line.strip()

    fastafile.close()

    return fastaSeqDict


def fixSpliceSense(fastaSequence, chrom, startRight, stopLeft, sense=""):
    spliceSense = {"GTAG": "+",
                   "GCAG": "+",
                   "ATAC": "+",
                   "CTAC": "-",
                   "CTGC": "-",
                   "GTAT": "-"
    }

    noncanonical = 0
    intronstart = stopLeft
    intronlen = startRight - stopLeft
    leftJunctionSig =fastaSequence[intronstart:intronstart+2]
    rightJunctionSig = fastaSequence[intronstart+intronlen-2:intronstart+intronlen]
    spliceJunction = leftJunctionSig + rightJunctionSig
    spliceJunction = spliceJunction.upper()
    if spliceSense.has_key(spliceJunction):
        sense = spliceSense[spliceJunction]
    else:
        noncanonical += 1
        senses = ["+", "-"]
        random.shuffle(senses)
        sense = senses[0]

    return sense, noncanonical


if __name__ == "__main__":
    main(sys.argv)