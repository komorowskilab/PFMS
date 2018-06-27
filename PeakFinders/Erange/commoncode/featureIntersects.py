#
#  featureIntersects.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
	pass

import sys
from cistematic.core import featuresIntersecting

print '%s: version 1.0' % sys.argv[0]

if len(sys.argv) < 5:
    print 'usage: python %s cistype tabfile [radius]' % sys.argv[0]
    sys.exit(1)

if len(sys.argv) > 1:
	cistype = sys.argv[1]
else:
	cistype = 'TFBSCONSSITES'

tabfile = open(sys.argv[2])

if len(sys.argv) > 3:
	radius = int(sys.argv[3])
else:
	radius = 100

previous = ''

posList = []
for line in tabfile:
        fields = line.split('\t')
        current = fields[0]
        if previous == current:
                continue
        previous = current
	chrom = fields[1][3:]
	posList.append((chrom, (int(fields[2]) + int(fields[3]))/2))

feats = featuresIntersecting('human', posList, radius, cistype)
featkeys = feats.keys()
featkeys.sort()
for (chrom, pos) in featkeys:
	print 'chr%s:%d-%d\t%s' % (chrom, pos, pos + 20, str(feats[(chrom, pos)]))
