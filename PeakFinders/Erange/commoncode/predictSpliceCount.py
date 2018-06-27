try:
	import psyco
	psyco.full()
except:
	print 'psyco not running'

from commoncode import *
import sys, math
print '%s: version 1.1' % sys.argv[0]

if len(sys.argv) < 6:
	print 'usage: python %s genome maxBorder uniquecountfile splicecountfile outfile' % sys.argv[0]
	sys.exit(1)

from cistematic.genomes import Genome
genome = sys.argv[1]
# number of nucleotides at the end of each exon that is affected by splicing
splicelead = int(sys.argv[2])
uniquefilecount = sys.argv[3]
splicefilecount =  sys.argv[4]
outfilename = sys.argv[5]

hg = Genome(genome)

gidDict = {}
gidList = []
uniqueCountDict = {}
spliceCountDict = {}

uniquefile = open(uniquefilecount)
for line in uniquefile:
	fields = line.strip().split()
	gidDict[fields[0]] = fields[1]
	gidList.append(fields[0])
	uniqueCountDict[fields[0]] = int(fields[2])
splicefile = open(splicefilecount)
for line in splicefile:
	fields = line.strip().split()
	spliceCountDict[fields[0]] = int(fields[2])

outfile = open(outfilename,'w')

gidList.sort()
for gid in gidList:
	symbol = gidDict[gid]
	featureList = hg.getGeneFeatures((genome, gid))
	newfeatureList = []
	featuresizesum = 0
	for (ftype, chrom, start, stop, sense) in featureList:
		if (start, stop) not in newfeatureList:
			newfeatureList.append((start, stop))
			featuresizesum += stop - start + 1
	if featuresizesum < 1:
		featuresizesum = 1
	splicearea = (len(newfeatureList) - 1) * splicelead
	if splicearea < splicelead:
		splicearea = 0
	fractionCoverage = featuresizesum / float(splicearea + featuresizesum)
	expectedSpliceCount = int(round(uniqueCountDict[gid]/fractionCoverage)) - uniqueCountDict[gid]
	# this p-value is based on the observed unique count, not the expected total count
	# nor the multi-read adjusted count
	pvalue = 1 - pow(1 - float(splicelead)/featuresizesum, uniqueCountDict[gid])
	print '%s %s %f %d %d' % (gid, symbol, pvalue, expectedSpliceCount, spliceCountDict[gid])
	outfile.write('%s\t%s\t%f\t%d\t%d\n' % (gid, symbol, pvalue, expectedSpliceCount, spliceCountDict[gid]))
#infile.close()
outfile.close()

