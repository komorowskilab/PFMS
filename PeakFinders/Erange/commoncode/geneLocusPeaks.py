#
#  geneLocusPeaks.py
#  ENRAGE
#

try:
	import psyco
	psyco.full()
except:
	pass

print '%s: version 2.0' % sys.argv[0]
from commoncode import *
import sys

if len(sys.argv) < 4:
	print 'usage: python %s genome rdsfile outfilename [-up upstream] [-down downstream] [-regions acceptfile] [-raw]' % sys.argv[0]
	print '\twhere upstream and downstream are in bp and and optional'
	sys.exit(1)

from cistematic.genomes import Genome
from cistematic.core.geneinfo import geneinfoDB

genome = sys.argv[1]
hitfile =  sys.argv[2]
outfilename = sys.argv[3]

upstream = 0
downstream = 0
useCDS = True
spanTSS = False
acceptfile = ''
acceptDict = {}
if '-up' in sys.argv:
    upstream = int(sys.argv[sys.argv.index('-up') + 1])

if '-down' in sys.argv:
    downstream = int(sys.argv[sys.argv.index('-down') + 1])

if '-regions' in sys.argv:
    acceptfile = sys.argv[sys.argv.index('-regions') + 1]
    acceptDict = getMergedRegions(acceptfile, maxDist=0, keepLabel=True, verbose=True)

if '-raw' in sys.argv:
    normalize = False
    normalizeBins = False
else:
    normalize = True
    normalizeBins = True    

doCache = False
if '-cache' in sys.argv:
    doCache = True

print 'upstream = %d downstream = %d' % (upstream, downstream)

hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
readlen = hitRDS.getReadSize()
normalizationFactor = 1.0
if normalize:
    totalCount = len(hitRDS)
    normalizationFactor = totalCount / 1000000.

hitDict = hitRDS.getReadsDict(doMulti=True, findallOptimize=True)

hg = Genome(genome)
idb = geneinfoDB(cache=True)

gidCount = {}
gidPos = {}
geneinfoDict = idb.getallGeneInfo(genome)
locusByChromDict = getLocusByChromDict(hg, upstream, downstream, useCDS=True, additionalRegionsDict=acceptDict)

gidList = hg.allGIDs()
gidList.sort()
for chrom in acceptDict:
    for (label, start, stop, length) in acceptDict[chrom]:
        if label not in gidList:
            gidList.append(label)
for gid in gidList:
    gidCount[gid] = 0

index = 0
for chrom in hitDict:
    if chrom not in locusByChromDict:
        continue
    print chrom
    for (start, stop, gid, glen) in locusByChromDict[chrom]:
        gidCount[gid] = 0.
        (topPos, numHits, smoothArray, numPlus) = findPeak(hitDict[chrom], start, glen, readlen)
        if len(topPos) > 0:
            gidCount[gid] = smoothArray[topPos[0]]
            gidPos[gid] = (chrom, start + topPos[0])
        else:
            gidPos[gid] = (chrom, start)
outfile = open(outfilename,'w')

for gid in gidList:
    if 'FAR' not in gid:
        symbol = 'LOC' + gid
        geneinfo = ''
        try:
            geneinfo = geneinfoDict[gid]
            symbol = geneinfo[0][0]
        except:
            pass
    else:
        symbol = gid
    if gid in gidCount and gid in gidPos:
        (chrom, pos) = gidPos[gid]
        outfile.write('%s\t%s\tchr%s\t%d\t%.2f\n' % (gid, symbol, chrom, pos, gidCount[gid]))
outfile.close()

