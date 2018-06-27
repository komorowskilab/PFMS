"""
MakeRdsFromBam

Created on Jun 3, 2010

@author: sau
"""

try:
    import psyco
    psyco.full()
except:
    pass

import sys, string, optparse, re
import pysam
from commoncode import readDataset, writeLog

verstring = "%prog: version 1.0"


def main(argv=None):
    if not argv:
        argv = sys.argv
    
    print verstring

    usage = "usage:  %prog label samfile outrdsfile [propertyName::propertyValue] [options]\
            \ninput reads must be sorted to properly record multireads"

    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--append", action="store_false", dest="init",
                      help="append to existing rds file [default: create new]")
    parser.add_option("--RNA", action="store_true", dest="rnaDataType",
                      help="set data type to RNA [default: DNA]")
    parser.add_option("-S", "--sam", action="store_true", dest="useSamFile",
                      help="input file is in sam format")
    parser.add_option("--index", action="store_true", dest="doIndex",
                      help="index the output rds file")
    parser.add_option("--cache", type="int", dest="cachePages",
                      help="number of cache pages to use [default: 100000")
    parser.add_option("-m", "--multiCount", type="int", dest="maxMultiReadCount",
                      help="multi counts over this value are discarded [default: 10]")
    parser.add_option("--rawreadID", action="store_false", dest="trimReadID",
                      help="use the raw read names")
    parser.set_defaults(init=True, doIndex=False, useSamFile=False, cachePages=100000,
                        maxMultiReadCount=10, rnaDataType=False, trimReadID=True)

    (options, args) = parser.parse_args(argv[1:])

    try:
        label = args[0]
    except IndexError:
        print "no label specified - see --help for usage"
        sys.exit(1)

    try:
        samFileName = args[1]
    except IndexError:
        print "no samfile specified - see --help for usage"
        sys.exit(1)

    try:
        outDbName = args[2]
    except IndexError:
        print "no outrdsfile specified - see --help for usage"
        sys.exit(1)

    makeRdsFromBam(label, samFileName, outDbName, options.init, options.doIndex, options.useSamFile,
                   options.cachePages, options.maxMultiReadCount, options.rnaDataType, options.trimReadID)


def makeRdsFromBam(label, samFileName, outDbName, init=True, doIndex=False, useSamFile=False,
                   cachePages=100000, maxMultiReadCount=10, rnaDataType=False, trimReadID=True):

    if useSamFile:
        fileMode = "r"
    else:
        fileMode = "rb"

    try:
        samfile = pysam.Samfile(samFileName, fileMode)
    except ValueError:
        print "samfile index not found"
        sys.exit(1)

    if rnaDataType:
        dataType = "RNA"
    else:
        dataType = "DNA"

    writeLog("%s.log" % outDbName, verstring, string.join(sys.argv[1:]))

    rds = readDataset(outDbName, init, dataType, verbose=True)
    if not init and doIndex:
        try:
            if rds.hasIndex():
                rds.dropIndex()
        except:
            pass

    if "sam_mapped" not in rds.getMetadata():
        rds.insertMetadata([("sam_mapped", "True")])

    defaultCacheSize = rds.getDefaultCacheSize()

    if cachePages > defaultCacheSize:
        if init:
            rds.setDBcache(cachePages, default=True)
        else:
            rds.setDBcache(cachePages)

    propertyList = []
    for arg in sys.argv:
        if "::" in arg:
            (pname, pvalue) = arg.strip().split("::")
            propertyList.append((pname, pvalue))

    if len(propertyList) > 0:
        rds.insertMetadata(propertyList)

    countReads = {"unmapped": 0,
                  "total": 0,
                  "unique": 0,
                  "multi": 0,
                  "multiDiscard": 0,
                  "splice": 0
    }

    readsize = 0
    insertSize = 100000

    uniqueInsertList = []
    multiInsertList = []
    spliceInsertList = []

    processedEntryDict = {}
    uniqueReadDict = {}
    multiReadDict = {}
    spliceReadDict = {}

    samFileIterator = samfile.fetch(until_eof=True)

    for read in samFileIterator:
        if read.is_unmapped:
            countReads["unmapped"] += 1
            continue

        if readsize == 0:
            take = (0, 2, 3) # CIGAR operation (M/match, D/del, N/ref_skip)
            readsize = sum([length for op,length in read.cigar if op in take])
            if init:
                rds.insertMetadata([("readsize", readsize)])

        #Build the read dictionaries
        try:
            readSequence = read.seq
        except KeyError:
            readSequence = ""

        pairReadSuffix = getPairedReadNumberSuffix(read)
        readName = "%s%s%s" % (read.qname, readSequence, pairReadSuffix)
        if trimReadID:
            rdsEntryName = "%s:%s:%d%s" % (label, read.qname, countReads["total"], pairReadSuffix)
        else:
            rdsEntryName = read.qname

        if processedEntryDict.has_key(readName):
            if isSpliceEntry(read.cigar):
                if spliceReadDict.has_key(readName):
                    del spliceReadDict[readName]
            else:
                if uniqueReadDict.has_key(readName):
                    del uniqueReadDict[readName]

                if multiReadDict.has_key(readName):
                    (read, priorCount, rdsEntryName) = multiReadDict[readName]
                    count = priorCount + 1
                    multiReadDict[readName] = (read, count, rdsEntryName)
                else:
                    multiReadDict[readName] = (read, 1, rdsEntryName)
        else:
            processedEntryDict[readName] = ""
            if isSpliceEntry(read.cigar):
                spliceReadDict[readName] = (read,rdsEntryName)
            else:
                uniqueReadDict[readName] = (read, rdsEntryName)

        if countReads["total"] % insertSize == 0:
            for entry in uniqueReadDict.keys():
                (readData, rdsEntryName) = uniqueReadDict[entry]
                chrom = samfile.getrname(readData.rname)
                uniqueInsertList.append(getRDSEntry(readData, rdsEntryName, chrom, readsize))
                countReads["unique"] += 1

            for entry in spliceReadDict.keys():
                (readData, rdsEntryName) = spliceReadDict[entry]
                chrom = samfile.getrname(readData.rname)
                spliceInsertList.append(getRDSSpliceEntry(readData, rdsEntryName, chrom, readsize))
                countReads["splice"] += 1

            for entry in multiReadDict.keys():
                (readData, count, rdsEntryName) = multiReadDict[entry]
                chrom = samfile.getrname(readData.rname)
                if count > maxMultiReadCount:
                    countReads["multiDiscard"] += 1
                else:
                    multiInsertList.append(getRDSEntry(readData, rdsEntryName, chrom, readsize, weight=count)) 
                    countReads["multi"] += 1

            rds.insertUniqs(uniqueInsertList)
            rds.insertMulti(multiInsertList)
            uniqueInsertList = []
            uniqueReadDict = {}
            multiInsertList = []
            multiReadDict = {}
            if dataType == "RNA":
                rds.insertSplices(spliceInsertList)
                spliceInsertList = []
                spliceReadDict = {}

            print ".",
            sys.stdout.flush()
            processedEntryDict = {}

        countReads["total"] += 1

    if len(uniqueReadDict.keys()) > 0:
        for entry in uniqueReadDict.keys():
            (readData, rdsEntryName) = uniqueReadDict[entry]
            chrom = samfile.getrname(readData.rname)
            uniqueInsertList.append(getRDSEntry(readData, rdsEntryName, chrom, readsize))
            countReads["unique"] += 1

        rds.insertUniqs(uniqueInsertList)

    if len(multiReadDict.keys()) > 0:
        for entry in multiReadDict.keys():
            (readData, count, rdsEntryName) = multiReadDict[entry]
            chrom = samfile.getrname(readData.rname)
            if count > maxMultiReadCount:
                countReads["multiDiscard"] += 1
            else:
                multiInsertList.append(getRDSEntry(readData, rdsEntryName, chrom, readsize, weight=count))
                countReads["multi"] += 1

        countReads["multi"] += len(multiInsertList)

    if len(spliceReadDict.keys()) > 0 and dataType == "RNA":
        for entry in spliceReadDict.keys():
            (readData, rdsEntryName) = spliceReadDict[entry]
            chrom = samfile.getrname(readData.rname)
            spliceInsertList.append(getRDSSpliceEntry(readData, rdsEntryName, chrom, readsize))
            countReads["splice"] += 1

        rds.insertSplices(spliceInsertList)

    countString = "\n%d unmapped reads discarded" % countReads["unmapped"]
    countString += "\t%d unique reads" % countReads["unique"]
    countString += "\t%d multi reads" % countReads["multi"]
    countString += "\t%d multi reads count > %d discarded" % (countReads["multiDiscard"], maxMultiReadCount)
    if dataType == "RNA":
        countString += "\t%d spliced reads" % countReads["splice"]

    print countString.replace("\t", "\n")

    writeLog("%s.log" % outDbName, verstring, countString)

    if doIndex:
        print "building index...."
        if cachePages > defaultCacheSize:
            rds.setDBcache(cachePages)
            rds.buildIndex(cachePages)
        else:
            rds.buildIndex(defaultCacheSize)


def getRDSEntry(alignedRead, readName, chrom, readSize, weight=1):
    start = int(alignedRead.pos)
    stop = int(start+readSize)
    sense = getReadSense(alignedRead.is_reverse)
    try:
        mismatchTag = alignedRead.opt("MD")
        mismatches = getMismatches(mismatchTag, alignedRead.seq, sense)
    except KeyError:
        mismatches = ""

    return (readName, chrom, start, stop, sense, 1.0/weight, '', mismatches)


def getRDSSpliceEntry(alignedRead, readName, chrom, readSize):
    (readName, chrom, start, stop, sense, weight, flag, mismatches) = getRDSEntry(alignedRead, readName, chrom, readSize)
    startL, startR, stopL, stopR = getSpliceBounds(start, readSize, alignedRead.cigar)
    
    return (readName, chrom, startL, stopL, startR, stopR, sense, 1.0, "", mismatches)


def getPairedReadNumberSuffix(read):
    readSuffix = ""
    if not isPairedRead(read):
        return ""

    if read.is_read1:
        readSuffix = "/1"
    elif read.is_read2:
        readSuffix = "/2"

    return readSuffix


def isPairedRead(read):
    return read.is_proper_pair and (read.is_read1 or read.is_read2)


def isSpliceEntry(cigarTupleList):
    isSplice = False
    for operation,length in cigarTupleList:
        if operation == 3:
            isSplice = True
            break

    return isSplice


def getReadSense(reverse):
    if reverse:
        sense = "-"
    else:
        sense = "+"

    return sense


def getMismatches(mismatchTag, querySequence="", sense="+", logErrors=False):
    output = []
    deletionMarker = "^"
    position = 0

    lengths = re.findall("\d+", mismatchTag)
    mismatchSequences = re.findall("\d+([ACGTN]|\\^[ACGTN]+)", mismatchTag)

    for mismatchEntry in range(len(mismatchSequences)):
        mismatch = mismatchSequences[mismatchEntry]
        position = position + int(lengths[mismatchEntry])
        if string.find(mismatch, deletionMarker) == 0:
            continue

        try:
            if querySequence:
                genomicNucleotide = querySequence[position]
            else:
                genomicNucleotide = "N"

            if sense == "-":
                mismatch = getComplementNucleotide(mismatch)
                genomicNucleotide  = getComplementNucleotide(genomicNucleotide)

            elandCompatiblePosition = int(position + 1)
            output.append("%s%d%s" % (mismatch, elandCompatiblePosition, genomicNucleotide))
            position += 1
        except IndexError:
            if logErrors:
                errorMessage = "getMismatch IndexError; tag: %s, seq: %s, pos: %d" % (mismatchTag, querySequence, position)
                writeLog("MakeRdsFromBamError.log", "1.0", errorMessage)

            return ""

    return string.join(output, ",")


def getComplementNucleotide(nucleotide):
    complement = {"A": "T",
                  "T": "A",
                  "C": "G",
                  "G": "C",
                  "N": "N"
    }

    return complement[nucleotide]


def getSpliceBounds(start, readsize, cigarTupleList):
    stopR = int(start + readsize)
    offset = 0

    for operation,length in cigarTupleList:
        if operation == 3:
            stopL = int(start + offset)
            startR = int(stopL + length)

            return start, startR, stopL, stopR
        else:
            offset += length


if __name__ == "__main__":
    main(sys.argv)