#
#  utrChanges.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    pass

from commoncode import *
from cistematic.genomes import Genome

import sys
print '%s: version 1.3' % sys.argv[0]

if len(sys.argv) < 4:
    print 'usage: python %s genome acceptedfile outfile' % sys.argv[0]
    sys.exit(1)

genome = sys.argv[1]
acceptfile =  sys.argv[2]
acceptDict = getMergedRegions(acceptfile, maxDist=0, keepLabel=True, verbose=True)
outfile = open(sys.argv[3],'w')

hg = Genome(genome)

origLocusByChromDict = getLocusByChromDict(hg, keepSense = True)
newLocusByChromDict = getLocusByChromDict(hg, additionalRegionsDict = acceptDict, keepSense = True)

new3utr = 0
new5utr = 0
changedGene = 0

for chrom in origLocusByChromDict:
    index = 0
    for (gstart, gstop, gid, glen, sense) in origLocusByChromDict[chrom]:   
        for (newstart, newstop, newgid, newlen, newsense) in newLocusByChromDict[chrom]:
            if gid == newgid:
                changedBoundary = False
                new3p = 'F'
                new5p = 'F'
                if newstart < gstart:
                    if sense == 'R':
                        new3utr += 1
                        new3p = 'T'
                        changedBoundary = True
                    elif sense == 'F':
                        new5utr += 1
                        new5p = 'T'
                        changedBoundary = True
                    else:
                        print sense
                if newstop > gstop:
                    if sense == 'R':
                        new5utr += 1
                        new5p = 'T'
                        changedBoundary = True
                    elif sense == 'F':
                        new3utr += 1
                        new3p = 'T'
                        changedBoundary = True
                    else:                
                        print sense
                if changedBoundary:
                    changedGene += 1
                    outfile.write('%s\tchr%s\t%d\t%d\t%s\tchr%s\t%d\t%d\t%s\t%s\n' % (gid, chrom, gstart, gstop, sense, chrom, newstart, newstop, new5p, new3p))
                continue
outfile.close()
print "%d new 5'utr" % new5utr
print "%d new 3'utr" % new3utr
print "%s affected genes" % changedGene

