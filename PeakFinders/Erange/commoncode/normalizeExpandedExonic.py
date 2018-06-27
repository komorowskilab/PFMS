try:
    import psyco
    psyco.full()
except:
    pass

import sys
print '%s: version 5.6' % sys.argv[0]
if len(sys.argv) < 6:
    print 'usage: python %s genome rdsfile uniqcountfile splicecountfile outfile [candidatefile acceptfile] [-gidField fieldID] [-maxLength kblength] [-cache]\n' % sys.argv[0]
    print "\twhere splicecountfile can be set to 'none' to not count splices\n"
    sys.exit(1)

from commoncode import *
from cistematic.genomes import Genome
from cistematic.core import chooseDB, cacheGeneDB, uncacheGeneDB

genome = sys.argv[1]
hitfile = sys.argv[2]
uniquecountfile = open(sys.argv[3])

dosplicecount = False
if sys.argv[4] != 'none':
    dosplicecount = True
    splicecountfile = open(sys.argv[4])

candidateLines = []
if len(sys.argv) > 6:
    try:
        candidatefile = open(sys.argv[6])
        candidateLines = candidatefile.readlines()
        candidatefile.close()
        acceptedfile = open(sys.argv[7],'w')
    except:
        pass

fieldID = 0
if '-gidField' in sys.argv:
    gidField = sys.argv.index('-gidField') + 1
    try:
        fieldID = int(sys.argv[gidField])
    except:
        print "couldn't process gidField"
        pass

maxLength = 1000000000.
if '-maxLength' in sys.argv:
    maxField = sys.argv.index('-maxLength') + 1
    try:
        maxLength = float(sys.argv[maxField])
    except:
        print "couldn't process maxLength"
        pass

extendGenome = ''
replaceModels = False
if '-models' in sys.argv:
    extendGenome = sys.argv[sys.argv.index('-models') + 1]
    if '-replacemodels' in sys.argv:
        replaceModels = True
        print "will replace gene models with %s" % extendGenome
    else:
        print "will extend gene models with %s" % extendGenome

doCache = False
if '-cache' in sys.argv:
    doCache = True
    cacheGeneDB(genome)
    hg = Genome(genome, dbFile=chooseDB(genome), inRAM=True)
    print '%s cached' % genome
else:
    hg = Genome(genome, inRAM=True)

if extendGenome != '':
    hg.extendFeatures(extendGenome, replace = replaceModels)
    
RDS = readDataset(hitfile, verbose = True, cache=doCache, reportCount=False)    
uniqcount = RDS.getUniqsCount()
print '%d unique reads' % uniqcount

splicecount = 0
countDict = {}
gidList = []
farList = []
candidateDict = {}


gidToGeneDict = {}
symbolToGidDict = {}
annotToGidDict = {}

featuresDict = hg.getallGeneFeatures()
print 'got featuresDict'

outfile = open(sys.argv[5],'w')

for line in uniquecountfile:
    fields = line.strip().split()
    gid = fields[fieldID]
    gene = fields[1]
    countDict[gid] = float(fields[-1])
    gidList.append(gid)
    gidToGeneDict[gid] = gene
uniquecountfile.close()

if dosplicecount:
    for line in splicecountfile:
        fields = line.strip().split()
        gid = fields[fieldID]
        try:
            countDict[gid] += float(fields[-1])
        except:
            print fields
            continue
        splicecount += float(fields[-1])
    splicecountfile.close()

for line in candidateLines:
    if '#' in line:
        continue
    fields = line.strip().split()
    gid = fields[1]
    gene = fields[0]
    if gid not in gidList:
        if gid not in farList:
            farList.append(gid)
            gidToGeneDict[gid] = gene
        if gid not in countDict:
            countDict[gid] = 0
        countDict[gid] += float(fields[6])
    if gid not in candidateDict:
        candidateDict[gid] = []
    candidateDict[gid].append((float(fields[6]), abs(int(fields[5]) - int(fields[4])), fields[3], fields[4], fields[5]))

totalCount = (uniqcount + splicecount) / 1000000.
uniqScale = uniqcount / 1000000.
for gid in gidList:
    gene = gidToGeneDict[gid]
    featureList = []
    try:
        featureList = featuresDict[gid]
    except:
        try:
            featureList = featuresDict[gene]
        except:
            print gene, gid
    newfeatureList = []
    geneLength = 0.
    for (ftype, chrom, start, stop, sense) in featureList:
        if (start, stop) not in newfeatureList:
            newfeatureList.append((start, stop))
            geneLength += (abs(start - stop) + 1.) / 1000.
    #rpkm is rpm / geneLength in kb.
    if geneLength < 0.1:
        geneLength = 0.1
    elif geneLength > maxLength:
        geneLength = maxLength
    rpm = countDict[gid] / totalCount
    rpkm = rpm / geneLength
    if gid in candidateDict:
        for (cCount, cLength, chrom, cStart, cStop) in candidateDict[gid]:
            cratio = cCount / (cLength / 1000.)
            cratio = (uniqScale * cratio) / totalCount
            if 10. * cratio < rpkm:
                continue
            countDict[gid] += cCount
            geneLength += cLength / 1000.
            acceptedfile.write('%s\t%s\t%s\t%s\t%.2f\t%d\t%s\n' % (gid, chrom, cStart, cStop, cratio, cLength, gene)) 
    rpm = countDict[gid] / totalCount
    rpmk = rpm / geneLength
    outfile.write('%s\t%s\t%.4f\t%.2f\n' %  (gid, gene, geneLength, rpkm))

for gid in farList:
    gene = gidToGeneDict[gid]
    geneLength = 0
    for (cCount, cLength, chrom, cStart, cStop) in candidateDict[gid]:
        geneLength += cLength / 1000.
    if geneLength < 0.1:
        continue
    for (cCount, cLength, chrom, cStart, cStop) in candidateDict[gid]:
        cratio = cCount / (cLength / 1000.)
        cratio = cratio / totalCount
        acceptedfile.write('%s\t%s\t%s\t%s\t%.2f\t%d\t%s\n' % (gene, chrom, cStart, cStop, cratio, cLength, gene)) 
    rpm = countDict[gid] / totalCount
    rpkm = rpm / geneLength
    outfile.write('%s\t%s\t%.4f\t%.2f\n' %  (gene, gene, geneLength, rpkm))
outfile.close()
try:
    acceptedfile.close()
except:
    pass

if '-cache' in sys.argv:
    uncacheGeneDB(genome)
