#
#  partition.py
#  ENRAGE
#
try:
    import psyco
    psyco.full()
except:
    pass

import sys, string
from commoncode import getMergedRegions, writeLog

versionString = '%s: version 2.0' % sys.argv[0]
print versionString

minFeature = 25
if len(sys.argv) < 4:
    print "usage: python %s mergeID regionfile1[,regionfile2,...] combpartitionfile [-minFeature bp] [-padregion bp] [-mergeregion bp] [-nomerge] [-log altlogfile] [-locid] [-ignorerandom] [-chromField fieldNum]\n"
    print "       where the regionfiles must be comma-separated with no white space"
    print "       -minFeature (default: %d bp) controls the size of the smallest partition\n" % minFeature
    sys.exit(1)

mergeID = sys.argv[1]
regionfiles = sys.argv[2]
outfilename = sys.argv[3]
logfilename = 'partition.log'
locID = False
ignoreRandom = False

if '-minFeature' in sys.argv:
    minFeature = int(sys.argv[sys.argv.index('-minFeature') + 1])

cField = 1
if '-chromField' in sys.argv:
        cField = int(sys.argv[sys.argv.index('-chromField') + 1])

padregion = 0
if '-padregion' in sys.argv:
    padField = sys.argv.index('-padregion') + 1
    padregion = int(sys.argv[padField])
    print "padding %d bp on each side of a region" % padregion

mergeregion = 0
if '-mergeregion' in sys.argv:
    mergeField = sys.argv.index('-mergeregion') + 1
    mergeregion = int(sys.argv[mergeField])
    print "merging regions closer than %d bp" % mergeregion

merging = True
if '-nomerge' in sys.argv:
    merging = False

if '-log' in sys.argv:
    logfilename = sys.argv[sys.argv.index('-log') + 1]

if '-locid' in sys.argv:
    locID = True
    print "using locations as region ID"

if '-norandom' in sys.argv:
    ignoreRandom = True
    print "ignoring 'random' chromosomes"

writeLog(logfilename, versionString, string.join(sys.argv[1:]))

allregionsDict = {}
regionFileList = regionfiles.split(',')
numRegions = len(regionFileList)
chromList = []
for regionID in range(numRegions):
    allregionsDict[regionID] = getMergedRegions(regionFileList[regionID], maxDist = mergeregion,  minHits=-1, fullChrom = True, verbose = True, chromField = cField, doMerge=merging, pad=padregion)
    for achrom in allregionsDict[regionID]:
        if achrom not in chromList:
            chromList.append(achrom)
            
outregionDict = {}

chromList = sorted(chromList)

for chrom in chromList:
    if ignoreRandom and 'random' in chrom:
        continue
    outregionDict[chrom] = []
    pointList = []
    for regionID in range(numRegions):
        if chrom in allregionsDict[regionID]:
            for (rstart, rstop, rlength) in allregionsDict[regionID][chrom]:
                pointList.append(rstart)
                pointList.append(rstop)
    pointList.sort()
    start = 0
    for point in pointList:
        if (point - start) > minFeature:
            outregionDict[chrom].append((start, point - 1, point - 1 - start))
            start = point

outfile = open(outfilename, 'w')
if locID:
    outfile.write('#chrom:start-stop\tchrom\tstart\tstop\tlength_kb\n')
else:
    outfile.write('#labelID\tchrom\tstart\tstop\tlength_kb\n')
index = 0
for chrom in outregionDict:
    for (start, stop, length) in outregionDict[chrom]:
        index += 1
        if locID:
            outfile.write("%s:%d-%d\t%s\t%d\t%d\t%.3f\n" % (chrom, start, stop, chrom, start, stop, length/1000.))
        else:
            outfile.write("%s%d\t%s\t%d\t%d\t%.3f\n" % (mergeID, index, chrom, start, stop, length/1000.))

message = "%s was partitioned into %d regions" % (mergeID, index)
print message
writeLog(logfilename, versionString, message)
