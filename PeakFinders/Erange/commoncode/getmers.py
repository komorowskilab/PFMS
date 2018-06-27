#
#  getmers.py
#  ENRAGE
#

import sys
try:
	import psyco
	psyco.full()
except:
	print 'psyco not running'

from cistematic.genomes import Genome

print '%s: version 1.1' % sys.argv[0]
if len(sys.argv) < 5:
	print 'usage: python %s genome merlen chrAny:start-stop outfile' % sys.argv[0]
	sys.exit(1)

genome = sys.argv[1]
merlen = int(sys.argv[2])
location = sys.argv[3]
outfilename = sys.argv[4]

(chrom, pos) = location.split(':')
chrom = chrom[3:]
(start, stop) = pos.split('-')
start = int(start)
regionlength = int(stop) - start + 1

hg = Genome(genome)

seq = hg.sequence(chrom, start, regionlength)

outfile = open(outfilename,'w')
print 'writing %d %d-mers' % (regionlength - merlen, merlen)
for index in range(regionlength - merlen):
    outfile.write(seq[index:index + merlen].upper() + '\n')

outfile.close()
