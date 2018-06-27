#
#  regionCounts.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    print 'psyco not running'

import sys
from commoncode import *

versionString = '%s: version 3.9' % sys.argv[0]
print versionString

if len(sys.argv) < 4:
    print 'usage: python %s regionfile rdsfile outfilename [-markRDS] [-chromField fieldNum] [-fullchrom] [-raw] [-padregion bp] [-mergeregion bp] [-nomerge] [-noUniqs] [-noMulti] [-splices] [-peak] [-cache [pages]] [-log altlogfile] [-rpkm] [-length] [-force]' % sys.argv[0]
    sys.exit(1)

regionfilename = sys.argv[1]
hitfile =  sys.argv[2]
outfilename = sys.argv[3]

logfilename = 'regionCounts.log'
if '-log' in sys.argv:
    logfilename = sys.argv[sys.argv.index('-log') + 1]

flagRDS = False
if '-markRDS' in sys.argv:
    flagRDS = True

cField = 1
if '-chromField' in sys.argv:
        cField = int(sys.argv[sys.argv.index('-chromField') + 1])

normalize = True
doRPKM = False
if '-rpkm' in sys.argv:
    doRPKM = True
elif '-raw' in sys.argv:
    normalize = False

doLength = False
if '-length' in sys.argv:
    doLength = True

useFullchrom = False
if '-fullchrom' in sys.argv:
    useFullchrom = True

forceRegion = False
if '-force' in sys.argv:
    forceRegion = True

padregion = 0
if '-padregion' in sys.argv:
    padField = sys.argv.index('-padregion') + 1
    padregion = int(sys.argv[padField])
    print 'padding %d bp on each side of a region' % padregion

mergeregion = 0
if '-mergeregion' in sys.argv:
    mergeField = sys.argv.index('-mergeregion') + 1
    mergeregion = int(sys.argv[mergeField])
    print 'merging regions closer than %d bp' % mergeregion

merging = True
if '-nomerge' in sys.argv:
    merging = False

usePeak = False
if '-peak' in sys.argv:
    print 'will use peak values'
    usePeak = True

doUniqs = True
doMulti = True
doSplices = False
if '-noUniqs' in sys.argv:
    doUniqs = False
if '-noMulti' in sys.argv:
    doMulti = False
if '-splices' in sys.argv:
    doSplices = True

doCache = False
cachePages = -1
if '-cache' in sys.argv:
    doCache = True
    try:
        cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])
    except: 
        pass

writeLog(logfilename, versionString, string.join(sys.argv[1:]))

regionDict = getMergedRegions(regionfilename, maxDist = mergeregion,  minHits=-1, keepLabel = True, fullChrom=useFullchrom, verbose = True, chromField = cField, doMerge=merging, pad=padregion)
labelList = []
labeltoRegionDict = {}
regionCount = {}

hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
readlen = hitRDS.getReadSize()
if cachePages > hitRDS.getDefaultCacheSize():
    hitRDS.setDBcache(cachePages)

totalCount = len(hitRDS)
if normalize:
    normalizationFactor = totalCount / 1000000.

chromList = hitRDS.getChromosomes(fullChrom=useFullchrom)
if len(chromList) == 0 and doSplices:
    chromList = hitRDS.getChromosomes(table='splices', fullChrom=useFullchrom)
chromList.sort()

if flagRDS:
    hitRDS.setSynchronousPragma('OFF')        

index = 0
for rchrom in regionDict:
    if forceRegion and rchrom not in chromList:
        print rchrom
        for (label, start, stop, length) in regionDict[rchrom]:
            regionCount[label] = 0
            labelList.append(label)
            labeltoRegionDict[label] = (rchrom, start, stop)

for rchrom in chromList:
    regionList = []
    if rchrom not in regionDict:
        continue
    print rchrom
    if useFullchrom:
        fullchrom = rchrom
    else:
        fullchrom = 'chr' + rchrom
    if usePeak:
        readDict = hitRDS.getReadsDict(chrom=fullchrom, withWeight=True, doMulti=True, findallOptimize=True)
        rindex = 0
        dictLen = len(readDict[fullchrom])
    
    for (label, start, stop, length) in regionDict[rchrom]:
        regionCount[label] = 0
        labelList.append(label)
        labeltoRegionDict[label] = (rchrom, start, stop)
        
    if useFullchrom:
        fullchrom = rchrom
    else:
        fullchrom = 'chr' + rchrom
    for (label, rstart, rstop, length) in regionDict[rchrom]:
        regionList.append((label, fullchrom, rstart, rstop))
        if usePeak:
            readList = []
            for localIndex in xrange(rindex, dictLen):
                read = readDict[fullchrom][localIndex]
                if read[0] < rstart:
                    rindex += 1
                elif rstart <= read[0] <= rstop:
                    readList.append(read)
                else:
                    break
            if len(readList) < 1:
                continue
            readList.sort()
            (topPos, numHits, smoothArray, numPlus) = findPeak(readList, rstart, rstop - rstart, readlen, doWeight=True)
            try:
                topValue = smoothArray[topPos[0]]
            except:
                print "problem with %s %s" % (str(topPos), str(smoothArray))
                continue
            regionCount[label] += topValue
        else:
            regionCount[label] += hitRDS.getCounts(fullchrom, rstart, rstop, uniqs=doUniqs, multi=doMulti, splices=doSplices)
    if flagRDS:
        hitRDS.flagReads(regionList, uniqs=doUniqs, multi=doMulti, splices=doSplices)
if flagRDS:
    hitRDS.setSynchronousPragma('ON')
            

#infile = open(infilename)
if normalize:
    for label in regionCount:
        regionCount[label] = float(regionCount[label]) / normalizationFactor

outfile = open(outfilename,'w')

if forceRegion:
	labelList.sort()
for label in labelList:
    (chrom, start, stop) = labeltoRegionDict[label]
    if useFullchrom:
        fullchrom = chrom
    else:
        fullchrom = 'chr' + chrom
    if normalize:
        if doRPKM:
            length = abs(stop - start) / 1000.
        else:
            length = 1.
        if length < 0.001:
            length = 0.001
        outfile.write('%s\t%s\t%d\t%d\t%.2f' % (label, fullchrom, start, stop, regionCount[label]/length))
        if doLength:
            outfile.write('\t%.1f' % length)
    else:
        outfile.write('%s\t%s\t%d\t%d\t%d' % (label, fullchrom, start, stop, regionCount[label]))
    outfile.write('\n')
#infile.close()
outfile.close()
if doCache and flagRDS:
    hitRDS.saveCacheDB(hitfile)

writeLog(logfilename, versionString, "returned %d region counts for %s (%.2f M reads)" % (len(labelList), hitfile, totalCount / 1000000.))
