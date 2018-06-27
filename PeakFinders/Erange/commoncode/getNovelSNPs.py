#
#  getNovelSNPs.py
#  ENRAGE
#
# This script attempts to annotate the novel sncs/snps from the snp summary file 
# Written by: Wendy Lee
# Written on: Aug 7th, 2008
# Last modified: Dec 14th, 2008 by Ali Mortazavi
import sys
from time import strftime
from cistematic.genomes import Genome 
from commoncode import *

print '%s: version 1.5' % sys.argv[0]

try:
    import psyco
    psyco.full()
except:
    pass

#Main program 
if len(sys.argv) < 3:
    print 'usage: python2.5 %s genome snpsfile nondbsnp_geneinfo_outfile' % sys.argv[0]
    sys.exit(1)

outStr = "" 
genome = sys.argv[1]
snpfile = sys.argv[2]
outfilename = sys.argv[3]

infile = file(snpfile, 'r')

hg = Genome(genome)
additionalDict={}

outS = ""
outfile  = open(outfilename, 'w')
outfile.write("#Sl\tCl\tchrom\tmis pos\t\tmatch\tuniq_mis\ttot_mis\tbase_chg\tknown_snp\tfunction\tgene\tgeneId\trpkm\n") 
for line in infile:
    if line[0] == '#':
        continue
    fields = line.split()
    if fields[8].find('N\A') == -1: 
        outfile.write(line)
    else:
        outS = ''
        gid = fields[11]
        snc_start = int(fields[3])
        featuresList = hg.getGeneFeatures((genome, gid))
        func = "N\A"
        for (ftype, chromosome, start, stop, orientation) in featuresList:
            if int(start) <= snc_start <= int(stop):
                func = ftype
                break  
        for i in range (0, 9):
            outS += fields[i] + "\t" 
        outS += func
        for i in range (10, 13):
            outS += "\t" + fields[i]
        outfile.write(outS +"\n")

outfile.close()

writeLog('snp.log', sys.argv[0], "outputfile: %s" % outfile)

