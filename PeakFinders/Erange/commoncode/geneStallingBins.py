#
#
#  geneStallingBins.py
#  ENRAGE
#

# originally from geneLocusBins.py
try:
    import psyco
    psyco.full()
except:
    pass

import sys
print '%s: version 1.3' % sys.argv[0]

from commoncode import *
from cistematic.genomes import Genome
from cistematic.core.geneinfo import geneinfoDB

if len(sys.argv) < 5:
    print 'usage: python %s genome rdsfile controlrdsfile outfilename [-upstream bp] [-downstream bp] [-regions acceptfile] [-cache] [-normalize] [-tagCount]' % sys.argv[0]
    sys.exit(1)

genome = sys.argv[1]
hitfile =  sys.argv[2]
controlfile = sys.argv[3]
outfilename = sys.argv[4]

normalizeBins = False
if '-normalize' in sys.argv:
    normalizeBins = True    

doTagCount = False
if '-tagCount' in sys.argv:
    doTagCount = True

acceptDict = {}
if '-regions' in sys.argv:
    acceptField = sys.argv.index('-regions') + 1
    acceptfile = sys.argv[acceptField]
    acceptDict = getMergedRegions(acceptfile, maxDist=0, keepLabel=True, verbose=True)

doCache = False
if '-cache' in sys.argv:
    doCache = True

bins = 4
if '-bins' in sys.argv:
    binfield = sys.argv.index('-bins') + 1
    bins = int(sys.argv[binfield])

upstreamBp = 300
downstreamBp = 0
doFlank = False

if '-upstream' in sys.argv:
    upfield = sys.argv.index('-upstream') + 1
    upstreamBp = int(sys.argv[upfield])
    doFlank = True

if '-downstream' in sys.argv:
    downfield = sys.argv.index('-downstream') + 1
    downstreamBp = int(sys.argv[downfield])
    doFlank = True

doCDS = True
limitNeighbor = False

hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
readlen = hitRDS.getReadSize()
hitNormalizationFactor = 1.0
if normalize:
    hitDictSize = len(hitRDS)
    hitNormalizationFactor = hitDictSize / 1000000.

controlRDS = readDataset(hitfile, verbose = True, cache=doCache)
controlNormalizationFactor = 1.0
if normalize:
    controlDictSize = len(hitRDS)
    controlNormalizationFactor = controlDictSize / 1000000.

hitDict = hitRDS.getReadsDict(doMulti=True, findallOptimize=True)
controlDict = controlRDS.getReadsDict(doMulti=True, findallOptimize=True)

hg = Genome(genome)
idb = geneinfoDB(cache=doCache)

gidBins = {}
gidLen = {}
geneinfoDict = idb.getallGeneInfo(genome)
locusByChromDict = getLocusByChromDict(hg, upstream=upstreamBp, downstream=downstreamBp, useCDS=doCDS, additionalRegionsDict=acceptDict, keepSense=True, adjustToNeighbor = limitNeighbor)

gidList = hg.allGIDs()
gidList.sort()
for chrom in acceptDict:
    for (label, start, stop, length) in acceptDict[chrom]:
        if label not in gidList:
            gidList.append(label)

(gidBins, gidLen) = computeRegionBins(locusByChromDict, hitDict, bins, readlen, gidList, hitNormalizationFactor, defaultRegionFormat=False, binLength=upstreamBp)
(controlBins, gidLen) = computeRegionBins(locusByChromDict, controlDict, bins, readlen, gidList, controlNormalizationFactor, defaultRegionFormat=False, binLength=upstreamBp)

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
        controlCount = 0.
        for binAmount in gidBins[gid]:
            tagCount += binAmount
        for binAmount in controlBins[gid]:
            controlCount += abs(binAmount)
        diffCount = tagCount + controlCount
        if diffCount < 0:
            diffCount = 0
	outfile.write('%s\t%s\t%.1f\t%d' % (gid, symbol, diffCount, gidLen[gid]))
        if (gidLen[gid] - 3 * upstreamBp) < upstreamBp:
                outfile.write('\tshort\n')
                continue
        TSSbins = (tagCount * (gidBins[gid][0] + gidBins[gid][1]) + controlCount * (controlBins[gid][0] + controlBins[gid][1])) / (upstreamBp / 50.)
        finalbin = (tagCount * gidBins[gid][-1] + controlCount * controlBins[gid][-1]) / ((gidLen[gid] - 3. * upstreamBp) / 100.)
        if finalbin <= 0.:
            finalbin = 0.01
        if TSSbins < 0:
            TSSbins = 0
        ratio =  float(TSSbins)/float(finalbin)
        for binAmount in gidBins[gid]:
            if doTagCount:
                binAmount = binAmount * tagCount / 100.
            if normalizeBins:
                if tagCount == 0:
                    tagCount = 1
                outfile.write('\t%.1f' % (100. * binAmount / tagCount))
            else:
                outfile.write('\t%.1f' % binAmount)
	outfile.write('\t%.2f\n' % ratio)
outfile.close()
