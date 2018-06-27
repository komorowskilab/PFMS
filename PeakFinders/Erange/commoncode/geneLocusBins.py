#
#  geneLocusBins.py
#  ENRAGE
#

# originally from version 1.3 of geneDownstreamBins.py
try:
    import psyco
    psyco.full()
except:
    pass

import sys
print '%s: version 2.1' % sys.argv[0]

from commoncode import *
from cistematic.genomes import Genome
from cistematic.core.geneinfo import geneinfoDB

if len(sys.argv) < 4:
    print 'usage: python %s genome rdsfile outfilename [-bins numbins] [-flank bp] [-upstream bp] [-downstream bp] [-nocds] [-regions acceptfile] [-cache] [-raw] [-force]' % sys.argv[0]
    sys.exit(1)

genome = sys.argv[1]
hitfile =  sys.argv[2]
outfilename = sys.argv[3]

normalizeBins = True
if '-raw' in sys.argv:
    normalizeBins = False    

acceptDict = {}
if '-regions' in sys.argv:
    acceptField = sys.argv.index('-regions') + 1
    acceptfile = sys.argv[acceptField]
    acceptDict = getMergedRegions(acceptfile, maxDist=0, keepLabel=True, verbose=True)

doCache = False
if '-cache' in sys.argv:
    doCache = True

bins = 10
if '-bins' in sys.argv:
    binfield = sys.argv.index('-bins') + 1
    bins = int(sys.argv[binfield])

upstreamBp = 0
downstreamBp = 0
doFlank = False
if '-flank' in sys.argv:
    flankfield = sys.argv.index('-flank') + 1
    upstreamBp = int(sys.argv[flankfield])
    downstreamBp = int(sys.argv[flankfield])
    doFlank = True

if '-upstream' in sys.argv:
    upfield = sys.argv.index('-upstream') + 1
    upstreamBp = int(sys.argv[upfield])
    doFlank = True

if '-downstream' in sys.argv:
    downfield = sys.argv.index('-downstream') + 1
    downstreamBp = int(sys.argv[downfield])
    doFlank = True

doCDS = True
if '-nocds' in sys.argv:
    doCDS = False
    
limitNeighbor = True
if '-force':
    limitNeighbor = False

hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
readlen = hitRDS.getReadSize()
normalizationFactor = 1.0
if normalizeBins:
    totalCount = len(hitRDS)
    normalizationFactor = totalCount / 1000000.

hitDict = hitRDS.getReadsDict(doMulti=True, findallOptimize=True)

hg = Genome(genome)
idb = geneinfoDB(cache=doCache)

gidBins = {}
gidLen = {}
geneinfoDict = idb.getallGeneInfo(genome)
if doFlank:
    locusByChromDict = getLocusByChromDict(hg, upstream=upstreamBp, downstream=downstreamBp, useCDS=doCDS, additionalRegionsDict=acceptDict, keepSense=True, adjustToNeighbor = limitNeighbor)
else:
    locusByChromDict = getLocusByChromDict(hg, additionalRegionsDict=acceptDict, keepSense=True)

gidList = hg.allGIDs()
gidList.sort()
for chrom in acceptDict:
    for (label, start, stop, length) in acceptDict[chrom]:
        if label not in gidList:
            gidList.append(label)

(gidBins, gidLen) = computeRegionBins(locusByChromDict, hitDict, bins, readlen, gidList, normalizationFactor, defaultRegionFormat=False)

outfile = open(outfilename,'w')

for gid in gidList:
    if 'FAR' not in gid:
        symbol = 'LOC' + gid
        geneinfo = ''
        try:
            geneinfo = geneinfoDict[gid]
            symbol = geneinfo[0][0]
        except:
            #print 'problem with %s' % (gid)
            pass
    else:
        symbol = gid
    if gid in gidBins and gid in gidLen:
        tagCount = 0.
        for binAmount in gidBins[gid]:
            tagCount += binAmount
	outfile.write('%s\t%s\t%.1f\t%d' % (gid, symbol, tagCount, gidLen[gid]))
	for binAmount in gidBins[gid]:
            if normalizeBins:
                if tagCount == 0:
                    tagCount = 1
                outfile.write('\t%.1f' % (100. * binAmount / tagCount))
            else:
                outfile.write('\t%.1f' % binAmount)
	outfile.write('\n')
outfile.close()
