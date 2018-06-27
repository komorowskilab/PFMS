#
#  geneDownstreamBins.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    pass

# originally from version 1.3 of geneDnaDownstreamCounts.py
from commoncode import *
from cistematic.genomes import Genome
from cistematic.core.geneinfo import geneinfoDB
import sys

print '%s: version 2.0' % sys.argv[0]
if len(sys.argv) < 4:
    print 'usage: %s genome rdsfile outfilename [-max regionSize]' % sys.argv[0]
    print '\n\twhere regionSize is the optional maximum region in bp\n'
    sys.exit(1)

genome = sys.argv[1]
hitfile =  sys.argv[2]
outfilename = sys.argv[3]

standardMinDist = 3000
if '-max' in sys.argv:
    standardMinDist = int(sys.argv[sys.argv.index('-max') + 1])

bins = 10
standardMinThresh = standardMinDist / bins

hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
readlen = hitRDS.getReadSize()
normalizationFactor = 1.0
if normalize:
    hitDictSize = len(hitRDS)
    normalizationFactor = hitDictSize / 1000000.

hitDict = hitRDS.getReadsDict(doMulti=True, findallOptimize=True)

hg = Genome(genome)
idb = geneinfoDB(cache=True)

gidDict = {}
geneinfoDict = idb.getallGeneInfo(genome)
featuresDict = hg.getallGeneFeatures()

#infile = open(infilename)
outfile = open(outfilename,'w')

gidList = hg.allGIDs()
gidList.sort()
for gid in gidList:
    symbol = 'LOC' + gid
    geneinfo = ''
    featureList = []
    try:
        geneinfo = geneinfoDict[gid]
        featureList = featuresDict[gid]
        symbol = geneinfo[0][0]
    except:
        print gid
    if len(featureList) == 0:
        continue
    newfeatureList = []
    for (ftype, chrom, start, stop, fsense) in featureList:
        if (start, stop) not in newfeatureList:
            newfeatureList.append((start, stop))
    if chrom not in hitDict:
        continue
    newfeatureList.sort()
    if len(newfeatureList) < 1:
        #print '%s %s %d' % (gid, symbol, -1)
        #outfile.write('%s\t%s\t%d\n' % (gid, symbol, -1))
        continue
    glen = standardMinDist
    if fsense == 'F':
        nextGene = hg.rightGeneDistance((genome, gid), glen * 2)
        if nextGene < glen * 2:
            glen = nextGene / 2
        if glen < 1:
            glen = 1
        gstart = newfeatureList[-1][1]
        gstop = newfeatureList[-1][1] + glen
    else:
        nextGene = hg.leftGeneDistance((genome, gid), glen * 2)
        if nextGene < glen * 2:
            glen = nextGene / 2
        if glen < 1:
            glen = 1
        gstart = newfeatureList[0][0] - glen
        if gstart < 0:
            gstart = 0
        gstop = newfeatureList[0][0]
    tagCount = 0
    if glen < standardMinDist:
        continue
    binList = [0.] * bins
    for (tagStart, sense, weight) in hitDict[chrom]:
        tagStart -= gstart
        if tagStart >= glen:
            break
        if tagStart > 0:
            tagCount += weight
            if fsense == 'F':
                # we are relying on python's integer division quirk
                binID = tagStart / standardMinThresh 
                binList[binID] += weight
            else:
                rdist = glen - tagStart
                binID = rdist / standardMinThresh 
                binList[binID] += weight
    if tagCount < 2:
        continue
    tagCount *= normalizationFactor
    print '%s %s %.2f %d %s' % (gid, symbol, tagCount, glen, str(binList))
    outfile.write('%s\t%s\t%.2f\t%d' % (gid, symbol, tagCount, glen))
    for binAmount in binList:
        outfile.write('\t%.2f' % binAmount)
    outfile.write('\n')
#infile.close()
outfile.close()

