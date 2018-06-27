#
#  regionBins.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    print 'psyco not running'

import sys
print '%s: version 2.0' % sys.argv[0]

if len(sys.argv) < 4:
    print 'usage: python %s regionfile rdsfile outfilename [-bins numbins] [-field fieldNum] [-raw] [-padregion bp] [-mergeregion bp] [-cache]' % sys.argv[0]
    sys.exit(1)

from commoncode import *

regionfilename = sys.argv[1]
hitfile =  sys.argv[2]
outfilename = sys.argv[3]

if '-raw' in sys.argv:
    normalize = False
    normalizeBins = False
else:
    normalize = True
    normalizeBins = True    

doCache = False
if '-cache' in sys.argv:
    doCache = True

cField = 1
if '-field' in sys.argv:
    fieldIndex = sys.argv.index('-field') + 1
    cField = int(sys.argv[fieldIndex])

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

bins = 10
if '-bins' in sys.argv:
    binfield = sys.argv.index('-bins') + 1
    bins = int(sys.argv[binfield])

hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
readlen = hitRDS.getReadSize()
normalizationFactor = 1.0
if normalize:
    totalCount = len(hitRDS)
    normalizationFactor = totalCount / 1000000.

chromList = hitRDS.getChromosomes(fullChrom=False)
chromList.sort()

regionDict = getMergedRegions(regionfilename, maxDist = mergeregion, keepLabel = True, verbose = True, chromField = cField, pad=padregion)

hitDict = hitRDS.getReadsDict(doMulti=True, findallOptimize=True)

(regionsBins, regionsLen) = computeRegionBins(regionsByChromDict, hitDict, bins, readlen, normalizationFactor)

outfile = open(outfilename, 'w')
for regionID in regionsBins:
    tagCount = 0.
    for binAmount in regionsBins[regionID]:
        tagCount += binAmount
    outfile.write('%s\t%s\t%.1f\t%d' % (regionID, regionID, tagCount, Len[gid]))
    for binAmount in gidBins[gid]:
            if normalizeBins:
                if tagCount == 0:
                    tagCount = 1
                outfile.write('\t%.1f' % (100. * binAmount / tagCount))
            else:
                outfile.write('\t%.1f' % binAmount)
    outfile.write('\n')

outfile.close()