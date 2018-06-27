#
#  makerdsfrombowtie.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 10/20/08.
#

try:
    import psyco
    psyco.full()
except:
    pass

import sys, string
from commoncode import readDataset, writeLog

verstring = "%s: version 4.1" % sys.argv[0]
print verstring

if len(sys.argv) < 4:
    print "usage: python %s label infilename outrdsfile [-RNA ucscGeneModels] [-append] [-index] [propertyName::propertyValue] [-spacer num] [-rawreadID] [-forcepair 1_or_2] [-flip] [-verbose] [-strip] [-cache numPages]" % sys.argv[0]
    sys.exit(1)

label = sys.argv[1]
filename = sys.argv[2]
outdbname = sys.argv[3]
stripSpace = False

delimiter = '|'
init = True
if '-append' in sys.argv:
    init = False

doIndex = False
if '-index' in sys.argv:
    doIndex = True

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

trimReadID = True
if '-rawreadID' in sys.argv:
    print "using raw readIDs"
    trimReadID = False

forcePair = False
forceID = 0
if '-forcepair' in sys.argv:
    forceID = int(sys.argv[sys.argv.index('-forcepair') + 1])
    forcePair = True

flip = False
if '-flip' in sys.argv:
    print "flipping read sense"
    flip = True

spacer = 2
if '-spacer' in sys.argv:
    spacer = int(sys.argv[sys.argv.index('-spacer') + 1])

if '-strip' in sys.argv:
    stripSpace = True

readsize = 0
maxBorder = 0
index = 0
insertSize = 100000

writeLog(outdbname + '.log', verstring, string.join(sys.argv[1:]))

def decodeMismatches(mString, rsense):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    
    output = []
    mismatches = mString.split(',')
    for mismatch in mismatches:
        (pos,change) = mismatch.split(':')
        (genNT, readNT) = change.split('>')
        if rsense == '-':
            readNT = complement[readNT]
            genNT  = complement[genNT]
        # for eland-compatibility, we are 1-based
        output.append('%s%d%s' % (readNT, int(pos)+1, genNT))
    
    return string.join(output, ',')

geneDict = {}
mapDict = {}
seenSpliceList = []
if dataType == 'RNA':
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
line = infile.readline()
if stripSpace:
    line = line.replace(' ','')
fields = line.split()
readsize = len(fields[5])
pairedTest = fields[0][-2:]
paired = False
if pairedTest in ['/1','/2'] or forcePair:
    print "assuming reads are paired"
    paired = True
readsizeString = str(readsize)
splicesizeString = readsizeString

print "read size: %d bp" % readsize
if init:
    rds.insertMetadata([('readsize', readsize)])
    if paired:
        rds.insertMetadata([('paired', 'True')])
if 'bowtie_mapped' not in rds.getMetadata():
    rds.insertMetadata([('bowtie_mapped', 'True')])

if dataType == 'RNA' and 'spacer' not in rds.getMetadata():
    rds.insertMetadata([('spacer', spacer)])
infile.close()

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
for line in infile:
    lIndex += 1
    if stripSpace:
        line = line.replace(' ','')
    fields = line.strip().split()
    readID = fields[0]
    if trimReadID:
        readID = string.join(readID.split(':')[1:],':')
    if readID != prevID:
        listlen = len(readList)
        if trimReadID:
            prevID = label + '-' + prevID
        if forcePair:
            prevID += '/%d' % forceID 
        if listlen == 1:
            (sense, chrom, start, mismatches) = readList[0]
            if flip:
                if sense == '+':
                    sense = '-'
                else:
                    sense = '+'
            if '|' not in chrom:
                stop = start + readsize
                uInsertList.append((prevID, chrom, start, stop, sense, 1.0, '', mismatches))
                uIndex += 1
            elif dataType == 'RNA':
                currentSplice = chrom
                (model, spliceID, regionStart) = currentSplice.split(delimiter)
                if model not in geneDict:
                    prevID = readID
                else:
                    (gsense, blockCount, transLength, chrom, chromstarts, blockSizes) = geneDict[model]
                    spliceID = int(spliceID)
                    spliceCount = 0
                    cumulative = 0
                    rstart = int(start) - spacer
                    lefthalf = maxBorder - rstart
                    if lefthalf < 1 or lefthalf > maxBorder:
                        prevID = readID
                    else:
                        righthalf = readsize - lefthalf
                        startL = int(regionStart)  + rstart
                        stopL = startL + lefthalf
                        startR = chromstarts[spliceID + 1]
                        stopR = chromstarts[spliceID + 1] + righthalf
                        sInsertList.append((prevID, chrom, startL, stopL, startR, stopR, sense, 1.0, '', mismatches))
                        sIndex += 1
        elif listlen > 1:
            prevID = prevID + '::' + str(listlen)
            mIndex += 1
            # ignore multireads that can also map across splices
            skip = False
            for (sense, chrom, start, mismatches) in readList:
                if '|' in chrom:
                    skip = True
            if not skip:
                for (sense, chrom, start, mismatches) in readList:
                    stop = start + readsize
                    if flip:
                        if sense == '+':
                            sense = '-'
                        else:
                            sense = '+'
                    mInsertList.append((prevID, chrom, start, stop, sense, 1.0 / listlen, '', mismatches))
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
        index += 1
    # add the new read
    sense = fields[1]
    chrom = fields[2]
    # for eland compat, we are 1-based
    start = int(fields[3]) + 1
    mismatches = ''
    if ':' in fields[-1]:
        mismatches = decodeMismatches(fields[-1], sense)
    readList.append((sense, chrom, start, mismatches))
    if lIndex % 1000000 == 0:
        print 'processed %d lines' % lIndex

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

