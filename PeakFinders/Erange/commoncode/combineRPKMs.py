#
#  combineRPKMS.py
#  ENRAGE
#

print 'version 1.0'
try:
	import psyco
	psyco.full()
except:
	pass

import sys

if len(sys.argv) < 2:
    print 'usage: python %s firstRPKM expandedRPKM finalRPKM combinedOutfile [-withmultifraction]' % sys.argv[0]
    sys.exit(1)

firstfile = open(sys.argv[1])
expandedfile = open(sys.argv[2])
finalfile = open(sys.argv[3])
outfile = open(sys.argv[4], 'w')

doFraction = False
if '-withmultifraction' in sys.argv:
    doFraction = True

firstDict = {}
gidDict = {}
expandedDict = {}

for line in firstfile:
    fields = line.strip().split()
    firstDict[fields[1]] = fields[-1]
firstfile.close()

for line in expandedfile:
    fields = line.strip().split()
    expandedDict[fields[1]] = fields[-1]
    gidDict[fields[1]] = fields[0]
expandedfile.close()

header ='gid\tRNAkb\tgene\tfirstRPKM\texpandedRPKM\tfinalRPKM'
if doFraction:
    header += '\tfractionMulti'
outfile.write(header + '\n')

for line in finalfile:
    fields = line.strip().split()
    gene = fields[0]
    rnakb = fields[1]
    finalRPKM = fields[2]
    if gene in firstDict:
        outline = '%s\t%s\t%s\t%s\t%s\t%s' % (gidDict[gene], rnakb, gene, firstDict[gene], expandedDict[gene], finalRPKM)
    else:
        outline = '%s\t%s\t%s\t%s\t%s\t%s' % (gidDict[gene], rnakb, gene, '', expandedDict[gene], finalRPKM)

    if doFraction:
        fraction = fields[3]
        outline += '\t' + fraction
    
    outfile.write(outline + '\n')

finalfile.close()
outfile.close()