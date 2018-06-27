try:
	import psyco
	psyco.full()
except:
	print 'psyco not running'

print 'version 3.6'

import sys

if len(sys.argv) < 3:
	print 'usage: python %s rdsfile outfilename [-cache pages]\n' % sys.argv[0]
	sys.exit(1)

from commoncode import *

hitfile =  sys.argv[1]
outfilename = sys.argv[2]
index = 0
startDict = {}
stopDict = {}
lengthDict = {}
resultDict = {}

doCache = False
cachePages = 100000
if '-cache' in sys.argv:
    doCache = True
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])

hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
if cachePages > hitRDS.getDefaultCacheSize():
    hitRDS.setDBcache(cachePages)

readlen = hitRDS.getReadSize()
hitDict = hitRDS.getSplicesDict(noSense=True)
outfile = open(outfilename,'w')

for chrom in hitDict:
    startDict[chrom] = []
    stopDict[chrom] = []
    resultDict[chrom] = []

index = 0
for chrom in hitDict:
    startFeature = 0
    for (tagStart, lstop, rstart, tagStop) in hitDict[chrom]:
        index += 1
        length = tagStop - tagStart
        if length < readlen + 5:
            continue
        startDict[chrom].append((tagStart, length))
        stopDict[chrom].append((tagStop, length))
    startDict[chrom].sort()
    stopDict[chrom].sort()

spliceEvent = 0
altSpliceEvent = 0
alternative = 1
for chrom in startDict:
    firstIndex = 0
    maxIndex = len(startDict[chrom])
    while firstIndex < maxIndex:
        (fstart, flen) = startDict[chrom][firstIndex]
        (start, length) = (fstart, flen)
        secondIndex = firstIndex
        secondLengths = []
        while (start - fstart) < readlen:
            if secondIndex >= maxIndex:
                break
            (start, length) = startDict[chrom][secondIndex]
            if (start - fstart) < readlen and abs(length - flen) > readlen:
                line =  (chrom, fstart, fstart + flen, chrom, start, start + length)
                alreadySeen = False
                for slength in secondLengths:
                    if abs(slength - length) < readlen:
                        alreadySeen = True
                if len(resultDict[chrom]) == 0:
                    resultDict[chrom].append(line)
                elif line != resultDict[chrom][-1] and not alreadySeen:
                    resultDict[chrom].append(line)
                    secondLengths.append(length)
                    altSpliceEvent += 1
                    spliceEvent += 1
            secondIndex += 1
        firstIndex = secondIndex
        spliceEvent += 1
    firstIndex = 0
    maxIndex = len(stopDict[chrom])
    while firstIndex < maxIndex:
        (fstop, flen) = stopDict[chrom][firstIndex]
        (stop, length) = (fstop, flen)
        secondIndex = firstIndex
        secondLengths = []
        while (stop - fstop) < readlen:
            if secondIndex >= maxIndex:
                break
            (stop, length) = stopDict[chrom][secondIndex]
            if (stop - fstop) < readlen and abs(length - flen) > readlen:
                line = (chrom, fstop - flen, fstop, chrom, stop - length, stop)
                alreadySeen = False
                for slength in secondLengths:
                    if abs(slength - length) < readlen:
                        alreadySeen = True
                if len(resultDict[chrom]) == 0:
                    resultDict[chrom].append(line)
                if line != resultDict[chrom][-1] and not alreadySeen:
                    resultDict[chrom].append(line)
                    secondLengths.append(length)
                    altSpliceEvent += 1
                    spliceEvent += 1
            secondIndex += 1
        firstIndex = secondIndex
        spliceEvent += 1
    resultDict[chrom].sort()
    for line in resultDict[chrom]:
        outfile.write('alt%d' % alternative + '\tchr%s\t%d\t%d\tchr%s\t%d\t%d\n'  % line)
        alternative += 1
    print chrom, maxIndex, spliceEvent, altSpliceEvent

print spliceEvent, altSpliceEvent
outfile.close()