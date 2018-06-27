#
#  listGeneFeatures.py
#  ENRAGE
#

import sys
print '%s: version 1.1' % sys.argv[0]
if len(sys.argv) < 4:
	print 'usage: python %s genome [acceptFile] gid outfile\n' % sys.argv[0]
	sys.exit(1)

from cistematic.genomes import Genome
from commoncode import *

genome = sys.argv[1]

additionalDict = {}
if len(sys.argv) == 4:
    gid = sys.argv[2]
    outfile = open(sys.argv[3],'w')
else:
    additionalDict = getMergedRegions(sys.argv[2], maxDist = 0, keepLabel = True, verbose = True)
    gid = sys.argv[3]
    outfile = open(sys.argv[4],'w')

hg = Genome(genome)

featuresDict = getFeaturesByChromDict(hg, additionalDict, restrictList=[gid])
outfile.write('track name="LOC%s"\n' % gid)

senseDict = {'F':'+', 'R':'-', '+':'+', '-':'-'}
for chrom in featuresDict:
    for (start, stop, fgid, sense, ftype) in featuresDict[chrom]:
	outfile.write('chr%s\t%d\t%d\t%s\t0\t%s\n' % (chrom, start, stop, ftype, senseDict[sense]))
outfile.close()

