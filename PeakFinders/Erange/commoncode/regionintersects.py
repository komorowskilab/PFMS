#
#  regionintersects.py
#  ENRAGE
#
try:
	import psyco
	psyco.full()
except:
	pass
from commoncode import *
import sys

print '%s: version 3.0' % sys.argv[0]

if len(sys.argv) < 6:
	print 'usage: python %s rdsfile1 regionfile1 rdsfile2 regionfile2 outfile [-reject1 File1] [-reject2 File2] [-union] [-cache] [-raw]' % sys.argv[0]
	sys.exit(1)

readOneName =  sys.argv[1]
regionOneName = sys.argv[2]
readTwoName = sys.argv[3]
regionTwoName = sys.argv[4]

#mergedist=200
mergedist=0

outfilename = sys.argv[5]
outfile = open(outfilename, 'w')

doCache = False
if '-cache' in sys.argv:
    doCache = True

trackReject = False
if '-union' in sys.argv:
    trackReject = True
    
doReject = False
if '-reject1' in sys.argv:
    trackReject = True
    doReject = True
    rejectOne = open(sys.argv[sys.argv.index('-reject1') + 1],'w')
if '-reject2' in sys.argv:
    trackReject = True
    doReject = True
    rejectTwo = open(sys.argv[sys.argv.index('-reject2') + 1],'w')

normalize = True
if '-raw' in sys.argv:
    normalize = False

doVerbose = False
if '-verbose' in sys.argv:
    doVerbose = True

oneDict = getMergedRegions(regionOneName, mergedist, verbose=doVerbose)
twoDict = getMergedRegions(regionTwoName, mergedist, verbose=doVerbose)

oneRDS = readDataset(readOneName, verbose=doVerbose, cache=doCache) 
twoRDS = readDataset(readTwoName, verbose=doVerbose, cache=doCache)

if normalize:
    normalize1 = len(oneRDS) / 1000000.
    normalize2 = len(twoRDS) / 1000000.
else:
    normalize1 = 1.
    normalize2 = 1.

commonRegions = 0
oneRejectIndex = 0
twoRejectIndex = 0

onePeaksDict = {}
oneFoundDict = {}
twoPeaksDict = {}

numRegionsOne = 0
numRegionsTwo = 0
for rchrom in oneDict:
    numRegionsOne += len(oneDict[rchrom])
for rchrom in twoDict:
    numRegionsTwo += len(twoDict[rchrom])
outfile.write('#%d\tregions in\t%s\n#%d\tregions in\t%s\n' % (numRegionsOne, regionOneName, numRegionsTwo, regionTwoName))

for rchrom in oneDict:
    if rchrom not in twoDict:
        continue
    print rchrom
    rindex = 0
    rindex2 = 0
    fullchrom = 'chr' + rchrom
    oneReads = oneRDS.getReadsDict(fullChrom=True, chrom=fullchrom, withWeight=True, doMulti=True)
    dictLen1 = len(oneReads[fullchrom])
    twoReads = twoRDS.getReadsDict(fullChrom=True, chrom=fullchrom, withWeight=True, doMulti=True)
    dictLen2 = len(twoReads[fullchrom])
    chrom = rchrom
    onePeaksDict[chrom] = []
    oneFoundDict[chrom] = []
    for (start, stop, length) in oneDict[chrom]:
        readList = []
        for localIndex in xrange(rindex, dictLen1):
            read = oneReads[fullchrom][localIndex]
            if read[0] < start:
                rindex += 1
            elif start <= read[0] <= stop:
                readList.append(read)
            else:
                break
        if len(readList) < 1:
            continue
        readList.sort()

        (topPos, numHits, smoothArray, numPlus) = findPeak(readList, start, length, doWeight=True)
        onePeakScore = smoothArray[topPos[0]]
        onePeaksDict[chrom].append((topPos[0] + start, length/2, start, stop, numHits/normalize1, onePeakScore/normalize1))
    for (start, stop, length) in twoDict[chrom]:
        readList2 = []
        for localIndex in xrange(rindex2, dictLen2):
            read = twoReads[fullchrom][localIndex]
            if read[0] < start:
                rindex2 += 1
            elif start <= read[0] <= stop:
                readList2.append(read)
            else:
                break
        if len(readList2) < 1:
            continue
        readList2.sort()
        
        (topPos, numHits, smoothArray, numPlus) = findPeak(readList2, start, length, doWeight=True)
        numHits /= normalize2
        twoIsCommon = False
        twoPeak = topPos[0] + start
        twoRadius = length/2
        twoPeakScore = smoothArray[topPos[0]] / normalize2
        for (onePeak, oneRadius, ostart, ostop, ohits, opeakScore) in onePeaksDict[chrom]:
            if abs(twoPeak - onePeak) < (twoRadius + oneRadius):
                if (onePeak, oneRadius, ostart, ostop, ohits) not in oneFoundDict:
                    oneFoundDict[chrom].append((onePeak, oneRadius, ostart, ostop, ohits))
                twoIsCommon = True
                commonRegions += 1
                outline = 'common%d\tchr%s\t%d\t%d\t%.1f\t%.1f\tchr%s\t%d\t%d\t%.1f\t%.1f' % (commonRegions, chrom, ostart, ostop, ohits, opeakScore, chrom, start, stop, numHits, twoPeakScore)
                if doVerbose:
                    print outline
                outfile.write(outline + '\n')
        if trackReject and not twoIsCommon:
            twoRejectIndex += 1
            outline = 'rejectTwo%d\tchr%s\t%d\t%d\t%.1f\t%.1f' % (twoRejectIndex, chrom, start, stop, numHits, twoPeakScore)
            if doReject:
                rejectTwo.write(outline + '\n')
            else:
                outfile.write(outline + '\n')
            if doVerbose:
                print outline

    if trackReject:
        #print len(onePeaksDict[chrom]), len(oneFoundDict[chrom])
        for (onePeak, oneRadius, ostart, ostop, ohits, opeakScore) in onePeaksDict[chrom]:
            if (onePeak, oneRadius, ostart, ostop, ohits) not in oneFoundDict[chrom]:
                oneRejectIndex += 1
                outline = 'rejectOne%d\tchr%s\t%d\t%d\t%.1f\t%.1f' % (oneRejectIndex, chrom, ostart, ostop, ohits, opeakScore)
                if doReject:
                    rejectOne.write(outline + '\n')
                else:
                    outfile.write(outline + '\n')
                if doVerbose:
                    print outline

if trackReject:
    print 'common: %d   one-only: %d   two-only: %d' % (commonRegions, oneRejectIndex, twoRejectIndex)
    outfile.write('#common: %d\tone-only: %d\ttwo-only: %d\n' % (commonRegions, oneRejectIndex, twoRejectIndex))
else:
    print 'common: %d' % commonRegions
    outfile.write('#common: %d\n' % commonRegions)

outfile.close()

