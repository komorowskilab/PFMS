#
#  makerdsfromeland2.py
#  ENRAGE
#
try:
    import psyco
    psyco.full()
except:
    pass

import sys, string
from commoncode import readDataset

print '%s: version 3.4' % sys.argv[0]
if len(sys.argv) < 4:
    print 'usage: python %s label infilename outrdsfile [-append] [-RNA ucscGeneModels] [propertyName::propertyValue] [-index] [-paired 1 or 2] [-extended] [-verbose] [-olddelimiter] [-maxlines num] [-cache numPages]' % sys.argv[0]
    sys.exit(1)

label = sys.argv[1]
filename = sys.argv[2]
outdbname = sys.argv[3]

delimiter = '|'
if '-olddelimiter' in sys.argv:
    delimiter = ':'

init = True
if '-append' in sys.argv:
    init = False

doIndex = False
if '-index' in sys.argv:
    doIndex = True

paired = False
pairID = '1'
if '-paired' in sys.argv:
    paired = True
    pairID = sys.argv[sys.argv.index('-paired') + 1]
    if pairID not in ['1','2']:
        print 'pairID value must be 1 or 2'
        sys.exit(-1)
    print 'Treating read IDs as paired with label = %s and pairID = %s' % (label, pairID)

dataType = 'DNA'
if '-RNA' in sys.argv:
    dataType = 'RNA'
    genedatafile = open(sys.argv[sys.argv.index('-RNA') + 1])

cachePages = 100000
if '-cache' in sys.argv:
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])

maxLines = 1000000000
if '-maxlines' in sys.argv:
    maxLines = int(sys.argv[sys.argv.index('-maxlines') + 1])

extended = False
if '-extended' in sys.argv:
    print 'using eland_extended input - will track mismatches'
    extended = True

verbose = False
if '-verbose' in sys.argv:
    verbose = True

readsize = 0
maxBorder = 0
index = 0
insertSize = 100000

def getUniqueMatch(elandCode):
    (zero, one, two) = elandCode.split(':')
    zero = int(zero)
    one = int(one)
    two = int(two)
    bestMatch = [False, False, False, False]
    if zero == 1:
        bestMatch[0] = True
        matchType = 0
    elif zero == 0 and one == 1:
        bestMatch[1] = True
        matchType = 1
    elif zero == 0 and one == 0 and two == 1:
        bestMatch[2] = True
        matchType = 2
    else:
        matchType = -1
    
    return (matchType, bestMatch)

def decodeMismatches(origSeq, code):
    output = []
    number = '0'
    index = 0
    for pos in code:
        if pos.isdigit():
            number += pos
        else:   
            index += int(number) + 1
            origNT = origSeq[index - 1]
            output.append('%s%d%s' % (origNT, index, pos))
            number = '0'
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
if cachePages > rds.getDefaultCacheSize():
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
        if pname == 'flowcell' and paired:
            pvalue = pvalue + '/' + pairID
        propertyList.append((pname, pvalue))
if len(propertyList) > 0:
    rds.insertMetadata(propertyList)

infile = open(filename,'r')
line = infile.readline()
fields = line.split()
readsize = len(fields[1])
readsizeString = str(readsize)
if dataType == 'RNA' and readsize > 32:
    splicesizeString = '32'
else:
    splicesizeString = readsizeString

print 'read size: %d bp' % readsize
if init:
    rds.insertMetadata([('readsize', readsize)])
    rds.insertMetadata([('eland_mapped', 'True')])
    if extended:
        rds.insertMetadata([('eland_extended', 'True')])
    if paired:
        rds.insertMetadata([('paired', 'True')])

trim = -4
if dataType == 'RNA':
    maxBorder = readsize + trim

insertList = []
infile = open(filename,'r')
print 'mapping unique reads...'
lineIndex = 0
for line in infile:
    lineIndex += 1
    if lineIndex > maxLines:
        break
    fields = line.split()
    if fields[2] in  ['QC','NM']:
        continue
    (matchType, bestMatch) = getUniqueMatch(fields[2])
    if matchType == -1:
        continue
    bestpos = []
    try:
        pos = fields[3].split(',')
    except:
        if verbose:
            print 'problem with line: %s' % line.strip()
        continue
    matchDict = {0:[], 1:[], 2:[], 3:[]}
    if len(pos) == 1:
        if 'splice' in pos:
            continue
        bestpos = pos
    else:
        currentChr = ''
        for apos in pos:
            if 'splice' in apos:
                continue
            if ':' in apos:
                (front, back) = apos.split(':')
                currentChr = front
            else:
                back = apos
                apos = currentChr + ':' + apos
            if extended:
                matchType = back.count('A') + back.count('C') + back.count('G') + back.count('T')
                if matchType > 2:
                    matchType = 3
            else:
                matchType = int(apos[-1])
            matchDict[matchType].append(apos)
            if bestMatch[matchType]:
                bestpos.append(apos)
    # for padded reads, mapped read might have more mismatches!
    if len(bestpos) == 0:
        # let's not worry about these yet.
        if 'splice' in line:
            continue
        for matchType in [1, 2, 3]:
            if len(matchDict[matchType]) > 0:
                if len(matchDict[matchType]) == 1 and 'splice' not in matchDict[matchType][0]:
                    bestpos = matchDict[matchType]
                break
        if len(bestpos) == 0 and verbose:
            print "couldn't pick best read from line: %s" % line
    
    for apos in bestpos:
        try:
            (chrom, back) = apos.split(':')
        except:
            continue
        if 'splice' in chrom:
            continue
        if '/' in chrom:
            chromfields = chrom.split('/')
            chrom = chromfields[-1]
        if '.' in chrom:
            try:
                (chrom, fileExt) = chrom.split('.')
            except:
                if verbose:
                    print 'problem with chromosome on line %s' % line.strip()
                continue
        if extended:
            if 'F' in back:
                sense = '+'
                (start, matchPart) = back.split('F')
            else:
                sense = '-'
                (start, matchPart) = back.split('R')
            start = int(start) 
            if matchPart == readsizeString:
                matchType = ''
            else:
                matchType = decodeMismatches(fields[1], matchPart)
        else:
            start = int(back[:-2])
            if back[-2] == 'F':
                sense = '+'		
            else:
                sense = '-'
        stop = int(start) + readsize - 1
        if paired:
            readID = label + '-' + str(lineIndex) + '/' + pairID
        else:
            readID = label + '-' + str(index)
        if len(chrom) > 0:
            insertList.append((readID, chrom, start, stop, sense, 1.0, '', matchType))
        if index % insertSize == 0:
            rds.insertUniqs(insertList)
            insertList = []
            print '.',
            sys.stdout.flush()
        index += 1

if len(insertList) > 0:
    rds.insertUniqs(insertList)
    insertList = []
print
print '%d unique reads' % index
infile.close()

if dataType == 'RNA':
    print 'mapping splices...'
    index = 0
    lineIndex = 0
    mapfile = open(filename,'r')
    for line in mapfile:
        lineIndex += 1
        if lineIndex > maxLines:
            break
        if 'splice' not in line:
            continue
        fields = line.strip().split()
        (matchType, bestMatch) = getUniqueMatch(fields[2])
        if matchType == -1:
            continue
        bestpos = []
        pos = fields[3].split(',')
        matchDict = {0:[], 1:[], 2:[], 3:[]}
        if len(pos) == 1:
            if 'chr' in pos:
                continue
            bestpos = pos
	else:
            currentSplice = ''
            for apos in pos:
                if 'splice' not in apos:
                    continue
                if ':' in apos:
                    if delimiter == ':':
                        try:
                            (extmodel, spliceID, regionStart, thepos) = apos.split(':')
                        except:
                            try:
                                (extmodel1, extmodel2, spliceID, regionStart, thepos) = apos.split(':')
                                extmodel = extmodel1 + ':' + extmodel2
                            except:
                                print 'warning: could not process splice %s' % apos
                                continue
                        currentSplice = extmodel + ':' + spliceID + ':' + regionStart
                    else:
                        try:
                            (currentSplice, thepos) = apos.split(':')
                        except:
                            try:
                                (extmodel1, restSplice, thepos) = apos.split(':')
                                currentSplice = extmodel1 + ':' + restSplice
                                (extmodel, spliceID, regionStart) = currentSplice.split(delimiter)
                            except:
                                print 'warning: could not process splice %s' % apos
                                continue
                else:
                    thepos = apos
                    apos = currentSplice + ':' + apos
                if extended:
                    matchType = thepos.count('A') + thepos.count('C') + thepos.count('G') + thepos.count('T')
                    if matchType > 2:
                        matchType = 3
                    # if readsize > 32, we risk loosing pefect matches that go beyond our expanded genome splices, so only ask for 32bp match
                    if thepos[:2] == splicesizeString:
                        matchType = 0
                else:
                    matchType = int(apos[-1])
                if bestMatch[matchType]:
                    bestpos.append(apos)
        # for padded reads, mapped read might have more mismatches!
        if len(bestpos) == 0:
            for matchType in [1, 2, 3]:
                if len(matchDict[matchType]) > 0:
                    if len(matchDict[matchType]) == 1 and 'splice' in matchDict[matchType][0]:
                        bestpos = matchDict[matchType]
                    break
            if len(bestpos) == 0 and verbose:
                print "couldn't pick best read from line: %s" % line
        for apos in bestpos:
            if delimiter == ':':
                try:
                    (extmodel, spliceID, regionStart, thepos) = apos.split(':')
                except:
                    try:
                        (extmodel1, extmodel2, spliceID, regionStart, thepos) = apos.split(':')
                        extmodel = extmodel1 + ':' + extmodel2
                    except:
                        print 'warning: could not process splice %s' % apos
                        continue
            else:
                try:
                    (currentSplice, thepos) = apos.split(':')
                except:
                    try:
                        (extmodel1, restSplice, thepos) = apos.split(':')
                        currentSplice = extmodel1 + ':' + restSplice
                    except:
                        print 'warning: could not process splice %s' % apos
                        continue
                (extmodel, spliceID, regionStart) = currentSplice.split(delimiter)
            modelfields = extmodel.split('/')
            if len(modelfields) > 2:
                model = string.join(modelfields[1:],'/')
            else:
                model = modelfields[1]
            if model not in geneDict:
                print fields
                continue
            (sense, blockCount, transLength, chrom, chromstarts, blockSizes) = geneDict[model]
            if extended:
                if 'F' in thepos:
                    rsense = '+'
                    (start, matchPart) = thepos.split('F')
                else:
                    rsense = '-'
                    (start, matchPart) = thepos.split('R')
                rstart = int(start) - 2 
                if matchPart == readsizeString:
                    matchType = ''
                elif matchPart[:2] == splicesizeString:
                    matchType = ''
                else:
                    matchType = decodeMismatches(fields[1], matchPart)
            else:
                rstart = int(thepos[:-2]) - 2
                if thepos[-2] == 'F':
                    rsense = '+'
                else:
                    rsense = '-'
            if trim <= rstart <= maxBorder:
                pass
            else:
                print rstart
                continue
            currentSplice = model + delimiter + spliceID + delimiter + regionStart
            spliceID = int(spliceID)
            spliceCount = 0
            cumulative = 0
            lefthalf = maxBorder - rstart
            if lefthalf < 1 or lefthalf > maxBorder:
                continue	
            righthalf = readsize - lefthalf
            startL = int(regionStart)  + rstart
            stopL = startL + lefthalf
            startR = chromstarts[spliceID + 1]
            stopR = chromstarts[spliceID + 1] + righthalf
            if paired:
                readName = label + '-' + str(lineIndex) + '/' + pairID
            else:
                readName = model + '-' + str(thepos)
            insertList.append((readName, chrom, startL, stopL, startR, stopR, rsense, 1.0, '', matchType))
            index += 1
            if index % insertSize == 0:
                rds.insertSplices(insertList)
                print '.',
                sys.stdout.flush()
                insertList = []
            if currentSplice not in seenSpliceList:
                seenSpliceList.append(currentSplice)
    mapfile.close()
    if len(insertList) > 0:
        rds.insertSplices(insertList)
        insertList = []
    print
    print 'saw %d spliced reads accross %d distinct splices' % (index, len(seenSpliceList))

infile = open(filename,'r')
print 'mapping multireads...'
lineIndex = 0
origReadid = rds.getMultiCount()
try:
    readid = int(origReadid) + 1
except:
    readid = 0
    origReadid = 0
print 'starting at %d' % (readid + 1)

for line in infile:
    lineIndex += 1
    if lineIndex > maxLines:
        break
    fields = line.split()
    if len(fields) < 4:
        continue
    if fields[2] == 'QC' or fields[2] == 'NM' or fields[3] == '-':
        continue
    (zero, one, two) = fields[2].split(':')
    zero = int(zero)
    one = int(one)
    two = int(two)
    
    bestMatch = [False] * readsize
    if zero > 1:
        bestMatch[0] = True
    elif zero == 0 and one > 1:
        bestMatch[1] = True
    elif zero == 0 and one == 0 and two > 1:
        bestMatch[2] = True
    else:
        continue
    #print '\n', lineIndex, line.strip()
    readcount = 0
    bestpos = []
    pos = fields[3].split(',')
    matchDict = {0:[], 1:[], 2:[], 3:[]}
    currentChr = ''
    for apos in pos:
        if ':' in apos:
            try:
                (front, back) = apos.split(':')
            except:
                if verbose:
                    print "problem splitting %s" % str(apos)
                continue
            currentChr = front
        else:
            back = apos
            apos = currentChr + ':' + apos
        if extended:
            matchType = back.count('A') + back.count('C') + back.count('G') + back.count('T')
        else:
            matchType = int(apos[-1])
        try:
            matchDict[matchType].append(apos)
        except:
            matchDict[matchType] = [apos]        
        if bestMatch[matchType]:
            bestpos.append(apos)
    # for padded reads, mapped read might have more mismatches!
    if len(bestpos) == 0:
        for matchType in [1, 2, 3]:
            if len(matchDict[matchType]) > 0:
                if len(matchDict[matchType]) > 1:
                    noSplice = True
                    for arg in matchDict[matchType]:
                        if 'splice' in arg:
                            noSplice = False
                    if noSplice:
                        bestpos = matchDict[matchType]
                break
        if len(bestpos) == 0 and verbose:
            print "couldn't pick best read from line: %s" % line
            continue
    hasSplice = False
    for apos in bestpos:
        if 'splice' in apos:
            hasSplice = True
    # do not allow multireads that can also map accross splices for now
    if hasSplice:
        if verbose:
            print "throwing out multiread because of splice conflict"
        continue
    
    if len(bestpos) > 0:
        readid += 1
        #print bestpos
    for apos in bestpos:
        readcount += 1
        (front, back) = apos.split(':')
        chrom = front[:-3]
        if extended:
            if 'F' in back:
                sense = '+'
                (start, matchPart) = back.split('F')
            else:
                sense = '-'
                (start, matchPart) = back.split('R')
            start = int(start)
            if matchPart == readsizeString:
                matchType = ''
            else:
                matchType = decodeMismatches(fields[1], matchPart)
        else:
            start = int(back[:-2])
            if back[-2] == 'F':
                sense = '+'
            else:
                sense = '-'
        stop = int(start) + readsize
        readName = '%dx%d' % (readid, len(bestpos))
        if paired:
            readName = label + '-' + str(lineIndex) + '/' + pairID + '::' + readName
        #print (readName, chrom, start, stop, sense, 1.0/len(bestpos), '', matchType)
        insertList.append((readName, chrom, start, stop, sense, 1.0/len(bestpos), '', matchType))
        if index % insertSize == 0:
            rds.insertMulti(insertList)
            insertList = []
            print '.',
            sys.stdout.flush()
        index += 1

if len(insertList) > 0:
    rds.insertMulti(insertList)
    insertList = []
print
print '%d multireads' % (readid - origReadid)

if doIndex:
    print 'building index....'
    rds.buildIndex(cachePages)

