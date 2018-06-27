#
#  transcripts.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 1/25/08.
#

import sys

tSize = 200000.
cellCount = 1e6
efficiency = 0.3

print '%s: version 3.0' % sys.argv[0]
if len(sys.argv) < 3:
    print 'usage: python %s rpkmFile outFile [-transcriptome size] [-cells count] [-efficiency fraction]' % sys.argv[0]
    print '\n\twhere transcriptome size is in Gbp (default %.1f), cell count is in arbitrary units (default %.1f) and efficiency is a fraction (default %.2f)\n' % (tSize, cellCount, efficiency) 
    sys.exit(1)

infile = open(sys.argv[1])
outfile = open(sys.argv[2],'w')

if '-transcriptome' in sys.argv:
    tSize = float(sys.argv[sys.argv.index('-transcriptome') + 1])

if '-cells' in sys.argv:
    cellCount = float(sys.argv[sys.argv.index('-cells') + 1])

if '-efficiency' in sys.argv:
    efficiency = float(sys.argv[sys.argv.index('-efficiency') + 1])

for line in infile:
    fields = line.strip().split()
    rpkm = float(fields[-1])
    transcripts = rpkm * tSize
    transPerCell = transcripts / cellCount / efficiency
    outfile.write('%s\t%.1f\t%.1f\n' % (fields[0], transcripts, transPerCell))
infile.close()
outfile.close()
