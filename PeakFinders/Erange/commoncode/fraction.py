#
#  fraction.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    pass

from random import random
import sys

print '%s: version 1.0' % sys.argv[0]

if len(sys.argv) < 4:
    print 'usage: python %s fraction infile outfile' % sys.argv[0]
    sys.exit(1)

fraction = float(sys.argv[1])
infile = open(sys.argv[2])
outfile = open(sys.argv[3], 'w')

totalIndex = 0
fractionIndex = 0
for line in infile:
    totalIndex += 1
    if random() <= fraction:
        outfile.write(line)
        fractionIndex += 1

infile.close()
outfile.close()

print '%d / %d = %.2f' % (fractionIndex, totalIndex, float(fractionIndex) / totalIndex)
