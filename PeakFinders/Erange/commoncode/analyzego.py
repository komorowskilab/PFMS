import sys
try:
	import psyco
	psyco.full()
except:
	print 'psyco not running'

import sys
print 'version 2.1'
if len(sys.argv) < 4:
    print 'usage: python %s genome infilename prefix [-geneName] [-field fieldID]' % sys.argv[0]
    sys.exit(1)

from cistematic.stat.analyzego import calculateGOStats
from cistematic.core.geneinfo import geneinfoDB

fieldID = 1
translateGene = False
if '-geneName' in sys.argv:
    translateGene = True
    fieldID = 0

if '-field' in sys.argv:
    fieldID = int(sys.argv[sys.argv.index('-field') + 1])

genome = sys.argv[1]
infilename = sys.argv[2]
prefix = sys.argv[3]

if translateGene:
    idb = geneinfoDB(cache=True)
    gidDict = {}
    geneinfoDict = idb.getallGeneInfo(genome)
    symbolToGidDict = {}
    for gid in geneinfoDict:
        symbol = geneinfoDict[gid][0][0].strip()
        symbolToGidDict[symbol] = gid

infile = open(infilename)

locusList = []
for line in infile:
    fields = line.split()
    if translateGene:
        gene = fields[fieldID]
        if 'LOC' in gene:
            gID = gene[3:]
        elif 'FAR' in gene:
            print 'ignoring %s' % gene
            continue
        else:
            try:
                gID = symbolToGidDict[gene]
            except:
                print 'ignoring %s' % gene
                continue
    else:
        gID = fields[fieldID]
    if (genome, gID) not in locusList:
        locusList.append((genome, gID))

if len(locusList) > 0:
    calculateGOStats(locusList, prefix)
