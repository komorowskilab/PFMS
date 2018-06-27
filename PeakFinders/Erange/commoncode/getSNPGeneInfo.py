#
#  getSNPGeneInfo.py
#  ENRAGE
#
# This script look for the gene info and expression level for the snps.
# Written by: Wendy Lee
# Written on: August 7th, 2008
# Last Modified:  December 14th, 2008 by Ali Mortazavi
try:
	import psyco
	psyco.full()
except:
	print 'psyco not running'
from cistematic.core import genesIntersecting, featuresIntersecting, cacheGeneDB, uncacheGeneDB
from cistematic.core.geneinfo import geneinfoDB
from cistematic.genomes import Genome
import sys

print '%s: version 4.5' % sys.argv[0]
if len(sys.argv) < 4:
    print 'usage: python %s genome snpsfile rpkmfile dbsnp_geneinfo_outfile [-cache] [-sense] [-flank bp]' % sys.argv[0]
    sys.exit(1)

nomatchfilename = ''
genome = sys.argv[1]
infilename = sys.argv[2]
rpkmfilename = sys.argv[3]
outfilename = sys.argv[4]
doCache = False
withSense = False
flankBP = 0
rpkmField = 3
    
if '-cache' in sys.argv:
    doCache = True
    cacheGeneDB(genome)
    idb = geneinfoDB(cache=True)
    print 'cached %s' % genome
else:
    idb = geneinfoDB()

if '-sense' in sys.argv:
    withSense = True

if '-flank' in sys.argv:
    flankBP = int(sys.argv[sys.argv.index('-flank')+1])

infile = open(infilename)
outfile = open(outfilename,'w')

geneinfoDict = idb.getallGeneInfo(genome)
posList = []
altPosDict = {}
altPosRevDict = {}
posLine = {}
altPosList = []

rpkmDict = {}
if rpkmfilename != 'NONE':
    rpkmfile = open(rpkmfilename, 'r')
    for l in rpkmfile:
        x = l.split()
        rpkmDict[x[0]] = x[rpkmField]
    rpkmfile.close()

index = 0
for line in infile:
    if line[0] == '#':
        continue
    fields = line.split('\t')
    chrom = fields[2][3:]
    start = int(fields[3])
    pos = (chrom, start)
    
    posList.append(pos)
    posLine[pos] = line
	
geneList = []
geneDict = {}
geneSense = {}

hg = Genome(genome)
#featureTypes = ['CDS'] + hg.getFeatureTypes('UT%')
featureTypes = ['CDS', 'UTR']
for ftype in featureTypes:
    if flankBP > 0:
        posDict = genesIntersecting(genome, posList, flank=flankBP)
    else:
        posDict = genesIntersecting(genome, posList)
    for pos in posDict:
        #print pos
        geneID = posDict[pos][0][0]
        try:
            symbol = geneinfoDict[geneID][0][0]
        except:
            symbol = 'LOC' + geneID
        try:
            geneDict[(symbol, geneID)].append(pos)
        except:
            geneList.append((symbol, geneID))
            geneDict[(symbol, geneID)] = [pos]
            geneSense[(symbol, geneID)] = posDict[pos][0][-1]
if doCache:
    uncacheGeneDB(genome)
outList = []
for (symbol, geneID) in geneList:
    geneDict[(symbol, geneID)].sort()
    seenLine = []
    for pos in geneDict[(symbol, geneID)]:
        if posLine[pos] in seenLine:
            continue
        try:
            rpkm = str(rpkmDict[geneID])
        except:
            rpkm = "N\A"
        if withSense:
            line = '%s\t%s\t%s\t%s\t%s' % (posLine[pos][:-1], symbol, geneID, rpkm, geneSense[(symbol,geneID)])
        else:
            line = '%s\t%s\t%s\t%s' % (posLine[pos][:-1], symbol, geneID, rpkm)
        seenLine.append(posLine[pos])
#        outfile.write(line+"\n")
        outList.append(line)

outList.sort(reverse=True)
for x in outList:
    outfile.write(x+"\n")
outfile.close()
