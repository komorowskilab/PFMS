#
#  makeSNPtrack.py
#  ENRAGE
#
# This script maps all the qualified SNC sites on to the genome browser 
# Output format: bed
# Written by: Wendy Lee
# Written on: August 18th, 2008
# Last Modified: December 14th, 2008 by Ali Mortazavi

import sys

print '%s: version 1.2' % sys.argv[0]

def mapSNPs(snpfile, track, outfilename):
    outfile = open(outfilename, 'w')
    baseColor = {'A':'200,0,255', 'T':'200,0,255', 'C':'200,0,255', 'G':'200,0,255'} 
    track ="track name=" + track + " description=" + track + " visibility=2 itemRgb=\"On\"\n"
    outfile.write(track)
    snpfile = open(snpfile, 'r')
    for line in snpfile:
        if line[0] == '#':
            continue
        fields = line.strip().split()
        chrom = fields[2]
        readstart = int(fields[3])-1
        readstop = readstart+1 
        base = fields[7]
        readName = base 
        color = baseColor[base[-1]] 
        if readName == 'A-G':
            color = '255,0,0'
        if readName == 'T-C':
            color = '0,0,255'
        score = str(0)
        sense = "+"
        outline = '%s\t%d\t%d\t%s\t%s\t%s\t-\t-\t\t%s\n' % (chrom, readstart, readstop, readName, score, sense, color)
        outfile.write(outline)
    snpfile.close()
    outfile.close()
	
#Main program 

if len(sys.argv) < 4:
    print 'usage: python2.5 %s snpfile trackname trackoutfile' % sys.argv[0]
    sys.exit(1)

snpfile = sys.argv[1]
track = sys.argv[2]
outfile = sys.argv[3]

mapSNPs(snpfile, track, outfile)
