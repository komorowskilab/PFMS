#
#  getfasta.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    pass

import sys
print "%s: version 3.4" % sys.argv[0]

from commoncode import *
from cistematic.genomes import Genome
from cistematic.core.geneinfo import geneinfoDB

if len(sys.argv) < 4:
    print "usage: python %s genome regionfile outfilename [-seqradius bp] [-minreads reads] [-returnTop num] [-maxsize bp] [-usepeak] [-dataset rdsfile] [-cache] [-compact]" % sys.argv[0]
    sys.exit(1)

genome = sys.argv[1]
regionfile = sys.argv[2]
outfilename = sys.argv[3]

# sequence radius
#seqsize = 500
seqsize = 50
minHitThresh = -1
topRegions = -1

if '-seqradius' in sys.argv:
    seqsize = int(sys.argv[sys.argv.index('-seqradius') + 1])

if '-minreads' in sys.argv:
    minHitThresh = int(sys.argv[sys.argv.index('-minreads') + 1])

doCache = False
if '-cache' in sys.argv:
    doCache = True

topRegions = 0
if '-returnTop' in sys.argv:
    topRegions = int(sys.argv[sys.argv.index('-returnTop') + 1])

maxsize = 300000000
if '-maxsize' in sys.argv:
    maxsize = int(sys.argv[sys.argv.index('-maxsize') + 1])
    
usePeaks = False
if '-usepeak' in sys.argv:
    usePeaks = True

doDataset = False
if '-dataset' in sys.argv:
    if usePeaks:
        print "ignoring dataset and relying on peak data"
    else:
        hitfile = sys.argv[sys.argv.index('-dataset') + 1]
        doDataset = True
        hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
        readlen = hitRDS.getReadSize()

doCompact = False
if '-compact' in sys.argv:
    doCompact = True

hg = Genome(genome)

outfile = open(outfilename,'w')

#readlen = readSize(hitfile)
#hitDict = getReadDict(hitfile)
if doCompact:
    regionDict = getMergedRegions(regionfile, minHits=minHitThresh, verbose=True, chromField = 0, compact = True, keepPeak=usePeaks, returnTop=topRegions)
else:
    regionDict = getMergedRegions(regionfile, minHits=minHitThresh, verbose=True, keepPeak=usePeaks, returnTop=topRegions)

ncregions = {}
for chrom in regionDict:
    ncregions[chrom] = []

index = 0


for achrom in regionDict:
    print "%s: processing %d regions" % (achrom, len(regionDict[achrom]))
    for region in regionDict[achrom]:
        border = 0
        if usePeaks:
            (rstart, rstop, rlen, peakPos, peakHeight) = region
            border = 200
        else:
            (rstart, rstop, rlen) = region
        if rlen > maxsize:
            print "%s:%d-%d length %d > %d max region size - skipping" % (achrom, rstart, rstop, rlen, maxsize)
            continue
        try:
            seq = hg.sequence(achrom, rstart - border, rlen + 2 * border)
        except:
            print "problem with %s" % str((rstart, rstop, rlen))
            continue
        if usePeaks:
            topPos = peakPos - rstart
            if peakHeight > minHitThresh:
                ncregions[achrom].append((rstart, rstop, rlen, [topPos], peakHeight))
                index += 1
        elif doDataset:
            thechrom = 'chr' + achrom
            print '.'
            hitDict = hitRDS.getReadsDict(chrom=thechrom, withWeight=True, doMulti=True, findallOptimize=True, start = rstart, stop = rstop)
            (topPos, numHits, smoothArray, numPlus) = findPeak(hitDict[thechrom], rstart, rlen, readlen)
            if numHits > minHitThresh:
                ncregions[achrom].append((rstart, rstop, rlen, topPos, numHits))
                index += 1
        else:
            ncregions[achrom].append((rstart, rstop, rlen, [-1], -1))
            index += 1
            

currentIndex = 1
for chrom in ncregions:
    for (rstart, rstop, rlen, topPos, numHits) in ncregions[chrom]:
        #print (chrom, rstart, rstop, rlen, str(topPos), numHits)
        if topPos[0] >= 0:
            newrstart = rstart + topPos[0] - seqsize
            newrlen = 2 * seqsize + 1
        else:
            newrstart = rstart
            newrlen = rlen
        seq2 = hg.sequence(chrom, newrstart, newrlen)
        #seq2 = hg.sequence(chrom, newrstart, 101, maskCDS=True, maskLower=True)
        outfile.write('>chr%s:%d-%d\n%s\n' % (chrom, newrstart, newrstart + newrlen, seq2))
outfile.close()
