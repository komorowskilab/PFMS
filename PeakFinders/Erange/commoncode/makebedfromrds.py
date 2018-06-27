#
#  makebedfromrds.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 7/19/08.
#

try:
    import psyco
    psyco.full()
except:
    pass

from commoncode import *
import sys

print '%s: version 3.0' % sys.argv[0]

plusColor = "0,0,255"
minusColor = "255,0,0"
multiPlusColor = "64,64,64"
multiMinusColor = "192,192,192"
spliceColor = "255,0,0"
uniqueColor = "0,0,0"
multiColor = "128,128,128"

if len(sys.argv) < 4:
    print 'usage: python %s trackLabel rdsfile bedfile [-nouniq] [-nomulti] [-splices] [-spliceColor] [-flag flagID] [-flagLike] [-pairs maxDist] [-cache pages] [-enforceChr] [-chrom chrID] [-strand + or -]' % sys.argv[0]
    sys.exit(1)

trackType = sys.argv[1]
rdsfile = sys.argv[2]
outfilename = sys.argv[3]

withUniqs = True
if '-nouniq' in sys.argv:
    withUniqs = False

withMulti = True
if '-nomulti' in sys.argv:
    withMulti = False

doSplices = False
if '-splices' in sys.argv:
    doSplices = True

doSpliceColor = False
if '-spliceColor' in sys.argv:
    doSpliceColor = True

doPairs = False
pairDist = 1000000
if '-pairs' in sys.argv:
    doPairs = True
    print 'will try to pair reads'
    try:
        pairDist = int(sys.argv[sys.argv.index('-pairs') + 1])
    except:
        pass

withFlag = ''
if '-flag' in sys.argv:
    withFlag = sys.argv[sys.argv.index('-flag') + 1]

useFlagLike = False
if '-flagLike' in sys.argv:
    useFlagLike = True

enforceChr = False
if '-enforceChr' in sys.argv:
    enforceChr = True

senseStrand = ''
if '-strand' in sys.argv:
    senseStrand = sys.argv[sys.argv.index('-strand') + 1]

allChrom = True
if '-chrom' in sys.argv:
    allChrom = False
    chromList = [sys.argv[sys.argv.index('-chrom') + 1]]

doCache = False
cachePages = 100000
if '-cache' in sys.argv:
    doCache = True
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])

if not withUniqs and not withMulti and not doSplices:
    print 'must be outputing at least one of uniqs, multi, or -splices - exiting'
    sys.exit(1)

print '\nsample:' 
RDS = readDataset(rdsfile, verbose = True, cache=doCache)

#check that this is better than the dataset's default cache size
if cachePages > RDS.getDefaultCacheSize():
    RDS.setDBcache(cachePages)

readlen = RDS.getReadSize()
minDist = -1 * readlen

if allChrom:
    if withUniqs:
        chromList = RDS.getChromosomes()
    elif withMulti:
        chromList = RDS.getChromosomes(table = 'multi')
    else:
        chromList = RDS.getChromosomes(table = 'splices')
    chromList.sort()

outfile = open(outfilename, 'w')
outfile.write('track name="%s" visibility=4 itemRgb="On"\n' % (trackType))

def singleReadWrite(chrom, pos, sense, weight, readID):
    start = pos
    stop = pos + readlen - 1
    if weight < 1.0:
        if sense == '+':
            senseCode = multiPlusColor
        else:
            senseCode = multiMinusColor
        outfile.write('%s %d %d %s %.1f %s 0 0 %s\n' % (chrom, start, stop, readID, weight, sense, senseCode))
    else:
        if sense == '+':
            senseCode = plusColor
        else:
            senseCode = minusColor
        outfile.write('%s %d %d %s %.1f %s 0 0 %s\n' % (chrom, start, stop, readID, weight, sense, senseCode))

def splitReadWrite(chrom, numPieces, startList, stopList, rsense, readName, plusSense, minusSense):
    readSizes = '%d' % (stopList[0] - startList[0])
    readCoords = '0'
    leftStart = startList[0]
    rightStop = stopList[-1]
    for index in range(1, numPieces):
        readSizes += ',%d' % (stopList[index] - startList[index])
        readCoords += ',%d' % (startList[index] - startList[0])
    if rsense == '+':
        senseCode = plusSense
    else:
        senseCode = minusSense
    outline = '%s\t%d\t%d\t%s\t1000\t%s\t0\t0\t%s\t%d\t%s\t%s\n' % (chrom, leftStart, rightStop, readName, rsense, senseCode, numPieces, readSizes, readCoords)
    outfile.write(outline)
    
for achrom in chromList:
    index = 0
    if achrom == 'chrM':
        continue
    if enforceChr and ('chr' not in achrom):
        continue
    print 'chromosome %s' % (achrom)
    # first handle uniqs and multireads
    if withUniqs or withMulti:
        if doPairs:
            hitDict = RDS.getReadsDict(fullChrom=True, chrom=achrom, flag=withFlag, withWeight=True, withPairID=True, doUniqs=withUniqs, doMulti=withMulti, readIDDict=True, flagLike=useFlagLike, strand=senseStrand)
            readIDList = hitDict.keys()
            if doSplices:
                spliceDict = RDS.getSplicesDict(fullChrom=True, chrom=achrom, flag=withFlag, withPairID=True, readIDDict=True, flagLike=useFlagLike, strand=senseStrand)
                spliceIDList = spliceDict.keys()
                combDict = {}
                for readID in readIDList:
                    combDict[readID] = 1
                for readID in spliceIDList:
                    combDict[readID] = 1
                combinedIDList = combDict.keys()
            else:
                combinedIDList = readIDList
            for readID in combinedIDList:
                localList = []
                try:
                    localList = hitDict[readID]
                except:
                    pass
                if doSplices:
                    try:
                        localList += spliceDict[readID]
                    except:
                        pass
                localList.sort()
                listLen = len(localList) - 1
                localIndex = 0
                while localIndex <= listLen:
                    try:
                        (leftpos, leftsense, leftweight, lPairID) = localList[localIndex]
                        leftstop = leftpos + readlen - 1
                        lpart = 1
                        startList = [leftpos]
                        stopList = [leftstop]
                    except:
                        (leftpos, LLstop, LRstart, leftstop, leftsense, lPairID) = localList[localIndex]
                        leftweight = 1.0
                        lpart = 2
                        startList = [leftpos, LRstart]
                        stopList = [LLstop, leftstop]
                    if localIndex < listLen:
                        try:
                            (rightpos, rightsense, rightweight, rPairID) = localList[localIndex + 1]
                            rightstop = rightpos + readlen - 1
                            rpart = 1
                            rstartList = [rightpos]
                            rstopList = [rightstop]
                        except:
                            (rightpos, RLstop, RRstart, rightstop, rightsense, rPairID) = localList[localIndex + 1]
                            rightweight = 1.0
                            rpart = 2
                            rstartList = [rightpos, RRstart]
                            rstopList = [RLstop, rightstop]
                    else:
                        rightsense = '+'
                        rightpos = 0
                        rstartList = []
                        rstopList = []
                    if leftsense == '+' and rightsense == '-' and minDist < (rightpos - leftstop) < pairDist and lPairID != rPairID:
                        if doSpliceColor:
                            if (lpart + rpart) > 2:
                                aColor = spliceColor
                                bColor = spliceColor
                            elif leftweight == 1.0 or rightweight == 1.0:
                                aColor = uniqueColor
                                bColor = uniqueColor
                            else:
                                aColor = multiColor
                                bColor = multiColor
                            splitReadWrite(achrom, lpart + rpart, startList + rstartList, stopList + rstopList, '+', readID, aColor, bColor)
                        elif leftweight == 1.0 or rightweight == 1.0:
                            splitReadWrite(achrom, lpart + rpart, startList + rstartList, stopList + rstopList, '+', readID, '0,0,0', minusColor)
                        else:
                            splitReadWrite(achrom, lpart + rpart, startList + rstartList, stopList + rstopList, '+', readID, '128,128,128', multiMinusColor)
                        localIndex += 2
                        index += 2
                    else:
                        if doSpliceColor:
                            if lpart  > 1:
                                aColor = spliceColor
                                bColor = spliceColor
                            elif leftweight == 1.0:
                                aColor = uniqueColor
                                bColor = uniqueColor
                            else:
                                aColor = multiColor
                                bColor = multiColor
                            splitReadWrite(achrom, lpart, startList, stopList, '+', readID, aColor, bColor)
                        elif leftweight == 1.0:
                            splitReadWrite(achrom, lpart, startList, stopList, leftsense, readID, plusColor, minusColor)
                        else:
                            splitReadWrite(achrom, lpart, startList, stopList, leftsense, readID, multiPlusColor, multiMinusColor)
                        localIndex += 1
                        index += 1
        else:
            hitDict = RDS.getReadsDict(fullChrom=True, chrom=achrom, flag=withFlag, withWeight=True, withID=True, doUniqs=withUniqs, doMulti=withMulti, readIDDict=False, flagLike=useFlagLike)
            try:
                for (pos, sense, weight, readID) in hitDict[achrom]:
                    #singleReadWrite(achrom, pos, sense, weight, readID)
                    splitReadWrite(achrom, 1, [pos], [pos + readlen - 1], sense, readID, plusColor, minusColor)
                    index += 1
            except:
                pass
            # now deal with splices
            if doSplices:
                spliceDict = RDS.getSplicesDict(fullChrom=True, chrom=achrom, flag=withFlag, withID=True, flagLike=useFlagLike)
                if achrom not in spliceDict:
                    continue
                for (readstart, Lstop, Rstart, readstop, rsense, readName) in spliceDict[achrom]:
                    splitReadWrite(achrom, 2, [readstart, Rstart], [Lstop, readstop], rsense, readName, plusColor, minusColor)
                    index += 1
    elif doSplices:
        spliceDict = RDS.getSplicesDict(fullChrom=True, chrom=achrom, flag=withFlag, withID=True, flagLike=useFlagLike)
        if achrom not in spliceDict:
            continue
        for (readstart, Lstop, Rstart, readstop, rsense, readName) in spliceDict[achrom]:
            splitReadWrite(achrom, 2, [readstart, Rstart], [Lstop, readstop], rsense, readName, plusColor, minusColor)
            index += 1
    print index
outfile.close()

    