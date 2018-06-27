#
#  geneLocusCounts.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    pass

import sys
print '%s: version 3.0' % sys.argv[0]
from commoncode import *

if len(sys.argv) < 4:
    print 'usage: python %s genome readDB outfilename [upstream] [downstream] [-noCDS] [-spanTSS] [-locusLength bplength] [-regions acceptfile] [-noUniqs] [-multi] [-splices]' % sys.argv[0]
    print '\twhere upstream and downstream are in bp and and optional'
    print '\tusing noCDS requires either upstream or downstream (but not both)'
    print '\tto be nonzero. Using -locuslength will report the first bplength'
    print '\tor the last bplength of the gene region depending on whether it'
    print '\tis positive or negative.\n'
    print '\twill by default only count the uniq reads (use -noUniqs to turn off)'
    print '\tbut can also count multi and splice reads given the appropriate flags\n'
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
if len(sys.argv) > 4:
    try:
        upstream = int(sys.argv[4])
    except:
        pass

if len(sys.argv) > 5:
    try:
        if '-' not in sys.argv[4]:
            downstream = int(sys.argv[5])
    except:
        pass

if '-noCDS' in sys.argv:
    useCDS = False

if '-spanTSS' in sys.argv:
    spanTSS = True

bplength = 0
if '-locusLength' in sys.argv:
    bplength = int(sys.argv[sys.argv.index('-locusLength') + 1])
    print 'returning only up to %d bp from gene locus' % bplength


if '-regions' in sys.argv:
    acceptfile = sys.argv[-1]
    acceptDict = getMergedRegions(acceptfile, maxDist=0, keepLabel=True, verbose=True)

print 'upstream = %d downstream = %d useCDS = %s spanTSS = %s' % (upstream, downstream, useCDS, spanTSS)
    
hitRDS = readDataset(hitfile, verbose = True)

readlen = hitRDS.getReadSize()

doUniqs = True
doMulti = False
doSplices = False
if '-noUniqs' in sys.argv:
    doUniqs = False
if '-multi' in sys.argv:
    doMulti = True
if '-splices' in sys.argv:
    doSplices = True

totalCount = hitRDS.getCounts(uniqs=doUniqs, multi=doMulti, splices=doSplices)

hg = Genome(genome)
idb = geneinfoDB(cache=True)

gidCount = {}
gidList = []
gidLen = {}
geneinfoDict = idb.getallGeneInfo(genome)
locusByChromDict = getLocusByChromDict(hg, upstream, downstream, useCDS, acceptDict, upstreamSpanTSS = spanTSS, lengthCDS = bplength)

locusChroms = locusByChromDict.keys()
chromList = hitRDS.getChromosomes(fullChrom=False)
chromList.sort()
for chrom in chromList:
    if chrom == 'M' or chrom not in locusChroms:
        continue
    print 'chr' + chrom
    fullchrom = 'chr' + chrom
    hitRDS.memSync(fullchrom, index=True)
    for (start, stop, gid, length) in locusByChromDict[chrom]:
        if gid not in gidList:
            gidList.append(gid)
            gidCount[gid] = 0
            gidLen[gid] = length
        gidCount[gid] += hitRDS.getCounts(fullchrom, start, stop, uniqs=doUniqs, multi=doMulti, splices=doSplices)
outfile = open(outfilename,'w')

totalCount /= 1000000.

outfile.write('#gid\tsymbol\tgidCount\tgidLen\trpm\trpkm\n')
gidList.sort()
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
    if gid in gidCount and gid in gidLen:
        rpm  = gidCount[gid] / totalCount
        rpkm = 1000. * rpm / gidLen[gid]
        outfile.write('%s\t%s\t%d\t%d\t%2.2f\t%2.2f\n' % (gid, symbol, gidCount[gid], gidLen[gid], rpm, rpkm))
outfile.close()

