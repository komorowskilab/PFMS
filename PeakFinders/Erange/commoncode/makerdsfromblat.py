#
#  makerdsfromblat.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 12/7/08.
#

try:
    import psyco
    psyco.full()
except:
    pass

import sys, string
from commoncode import readDataset, writeLog

verstring = "%s: version 3.9" % sys.argv[0]
print verstring

if len(sys.argv) < 4:
    print "usage: python %s label infilename outrdsfile [-append] [-index] [propertyName::propertyValue] [-rawreadID] [-forceRNA] [-flag] [-strict minSpliceLen] [-spliceonly] [-verbose] [-cache numPages]" % sys.argv[0]
    sys.exit(1)

label = sys.argv[1]
filename = sys.argv[2]
outdbname = sys.argv[3]

delimiter = '|'
init = True
if '-append' in sys.argv:
    init = False

doIndex = False
if '-index' in sys.argv:
    doIndex = True

minSpliceLength = 0
if '-strict' in sys.argv:
    minSpliceLength = int(sys.argv[sys.argv.index('-strict') + 1])
    print "requiring at least %d bp on each side of a splice" % minSpliceLength
minIntron = 10

dataType = 'DNA'
if '-RNA' in sys.argv:
    dataType = 'RNA'
    genedatafile = open(sys.argv[sys.argv.index('-RNA') + 1])

cachePages = 100000
if '-cache' in sys.argv:
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])

verbose = False
if '-verbose' in sys.argv:
    verbose = True

forceRNA = False
if '-forceRNA' in sys.argv:
    print "forcing datatype to RNA"
    dataType = 'RNA'
    forceRNA = True

theFlag = ''
if '-flag' in sys.argv:
    theFlag = 'blat'

spliceOnly = False
if '-spliceonly' in sys.argv:
    spliceOnly = True
    print "only recording splices"
    
trimReadID = True
if '-rawreadID' in sys.argv:
    print "using raw readIDs"
    trimReadID = False

readsize = 0
maxBorder = 0
index = 0
insertSize = 100000

writeLog(outdbname + '.log', verstring, string.join(sys.argv[1:]))

def decodeMismatches(gString, rString, rsense):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    
    output = []
    rlen = len(gString)
    partIndex = 0
    for rindex in xrange(rlen):
        if gString == ',':
            partIndex += 1
        if gString[rindex] == rString[rindex]:
            continue
        genNT = gString[rindex]
        readNT = rString[rindex]
        # for eland-compatibility, we are 1-based
        output.append('%s%d%s' % (readNT, rindex + 1 - partIndex, genNT))
            
    return string.join(output, ',')

geneDict = {}
mapDict = {}
seenSpliceList = []
if dataType == 'RNA' and not forceRNA:
    for line in genedatafile:
        fields = line.strip().split('\t')
        blockCount = int(fields[7])
        if blockCount < 2:
            continue
        uname = fields[0]
        chrom = fields[1]
        sense = fields[2]
        chromstarts = fields[8][:-1].split(',')
        chromstops = fields[9][:-1].split(',')
        exonLengths = []
        totalLength = 0
        for index in range(blockCount):
            chromstarts[index] = int(chromstarts[index])
            chromstops[index] = int(chromstops[index])
            exonLengths.append(chromstops[index] - chromstarts[index])
            totalLength += exonLengths[index]
        geneDict[uname] = (sense, blockCount, totalLength, chrom, chromstarts, exonLengths)
        mapDict[uname] = []
    genedatafile.close()

rds = readDataset(outdbname, init, dataType, verbose=True)

#check that our cacheSize is better than the dataset's default cache size
defaultCacheSize = rds.getDefaultCacheSize()
if cachePages > defaultCacheSize:
    if init:
        rds.setDBcache(cachePages, default=True)
    else:
        rds.setDBcache(cachePages)

if not init and doIndex:
    try:
        if rds.hasIndex():
            rds.dropIndex()
    except:
        if verbose:
            print "couldn't drop Index"

propertyList = []
for arg in sys.argv:
    if '::' in arg:
        (pname, pvalue) = arg.strip().split('::')
        propertyList.append((pname, pvalue))
if len(propertyList) > 0:
    rds.insertMetadata(propertyList)

# make some assumptions based on first read
infile = open(filename,'r')
for arg in range(6):
    line = infile.readline()
fields = line.split()
readsize = int(fields[10])
pairedTest = fields[9][-2:]
paired = False
if pairedTest in ['/1','/2']:
    print "assuming reads are paired"
    paired = True
readsizeString = fields[10]
splicesizeString = readsizeString

print "read size: %d bp" % readsize
if init:
    rds.insertMetadata([('readsize', readsize)])
    if paired:
        rds.insertMetadata([('paired', 'True')])
infile.close()
if 'blat_mapped' not in rds.getMetadata():
    rds.insertMetadata([('blat_mapped', 'True')])

minReadScore = readsize - readsize/25 - 1
trim = -4
if dataType == 'RNA':
    maxBorder = readsize + trim

infile = open(filename,'r')
prevID = ''
readList = []
uInsertList = []
mInsertList = []
sInsertList = []
index = uIndex = mIndex = sIndex = lIndex = 0
bestScore = 0
# skip headers
for arg in range(5):
    line = infile.readline()
for line in infile:
    lIndex += 1
    fields = line.strip().split()
    readID = fields[9]
    if trimReadID:
        readID = string.join(readID.split(':')[1:],':')
    if readID != prevID:
        newReadList = []
        if bestScore > minReadScore:
            for readData in readList:
                if readData[1] == bestScore:
                    newReadList.append(readData)
        if trimReadID:
            prevID = label + '-' + prevID
        listlen = len(newReadList)
        if listlen == 1:
            parts = int(newReadList[0][0])
            if parts == 1 and not spliceOnly:
                (part, score, sense, chrom, start, mismatches) = newReadList[0]
                stop = start + readsize
                uInsertList.append((prevID, chrom, start, stop, sense, 1.0, theFlag, mismatches))
                uIndex += 1
            elif forceRNA and parts == 2:
                (part, score, sense, chrom, startList, lengthList, mismatchList) = newReadList[0]
                startL = int(startList[0]) 
                stopL = startL + int(lengthList[0])
                startR = int(startList[1])
                stopR = startR + int(lengthList[1])
                if stopL + minIntron < startR:
                    sInsertList.append((prevID, chrom, startL, stopL, startR, stopR, sense, 1.0, theFlag, mismatches))
                    sIndex += 1
            elif parts == 2:
                print newReadList
                (part, score, sense, chrom, start, mismatches) = newReadList[0]
                currentSplice = chrom
                (model, spliceID, regionStart) = currentSplice.split(delimiter)
                if model not in geneDict:
                    print fields
                    continue
                (gsense, blockCount, transLength, chrom, chromstarts, blockSizes) = geneDict[model]
                spliceID = int(spliceID)
                spliceCount = 0
                cumulative = 0
                rstart = int(start) - 2
                lefthalf = maxBorder - rstart
                if lefthalf < 1 or lefthalf > maxBorder:
                    continue	
                righthalf = readsize - lefthalf
                startL = int(regionStart)  + rstart
                stopL = startL + lefthalf
                startR = chromstarts[spliceID + 1]
                stopR = chromstarts[spliceID + 1] + righthalf
                if stopL + minIntron < startR:
                    sInsertList.append((prevID, chrom, startL, stopL, startR, stopR, sense, 1.0, theFlag, mismatches))
                    sIndex += 1
        elif listlen > 1 and not spliceOnly:
            prevID = prevID + '::' + str(listlen)
            mIndex += 1
            # ignore multireads that can also map across splices
            skip = False
            for readData in newReadList:
                if readData[0] > 1:
                    skip = True
            if not skip:
                for (part, score, sense, chrom, start, mismatches) in newReadList:
                    stop = start + readsize
                    mInsertList.append((prevID, chrom, start, stop, sense, 1.0 / listlen, theFlag, mismatches))
        else:
            prevID = readID
        if index % insertSize == 0:
            rds.insertUniqs(uInsertList)
            rds.insertMulti(mInsertList)
            uInsertList = []
            mInsertList = []
            if dataType == 'RNA':
                rds.insertSplices(sInsertList)
                sInsertList = []
            print '.',
            sys.stdout.flush()
        # start processing new read
        readList = []
        prevID = readID
        bestScore = 0
        index += 1
    # add the new read
    score = int(fields[0])
    sense = fields[8]
    chrom = fields[13]
    parts = int(fields[17])
    passStrict = True
    if parts > 1:
        lengthList = fields[18][:-1].split(',')
        startList = fields[20][:-1].split(',')
        listlen = len(lengthList)
        for lpos in range(listlen):
            if int(lengthList[lpos]) < minSpliceLength:
                passStrict = False
            # throw out deletions, for now
            if lpos > 0:
                if int(lengthList[lpos - 1]) == int(startList[lpos]):
                    passStrict = False
        pass
    else:
        start = int(fields[15])
    if passStrict:
        if score > bestScore:
            bestScore = score
        mismatches = ''
        if int(fields[1]) > 0:
            try:
                mismatches = decodeMismatches(fields[-1].upper(), fields[-2].upper(), sense)
            except:
                mismatches = ''
        if parts == 1:
            readList.append((parts, score, sense, chrom, start, mismatches))
        else:
            readList.append((parts, score, sense, chrom, startList, lengthList, mismatches))
    if lIndex % 1000000 == 0:
        print "processed %d lines" % lIndex

print "%d lines processed" % lIndex
lines = []

if len(uInsertList) > 0:
    rds.insertUniqs(uInsertList)
if len(mInsertList) > 0:
    rds.insertMulti(mInsertList)
if len(sInsertList) > 0:
    rds.insertSplices(sInsertList)

combString = "%d unique reads" % uIndex
combString += "\t%d multi reads" % mIndex
if dataType == 'RNA':
    combString += "\t%d spliced reads" % sIndex

print
print combString.replace('\t', '\n')

writeLog(outdbname + '.log', verstring, combString)

if doIndex:
    print "building index...."
    if cachePages > defaultCacheSize:
        rds.setDBcache(cachePages)
        rds.buildIndex(cachePages)
    else:
        rds.buildIndex(defaultCacheSize)

