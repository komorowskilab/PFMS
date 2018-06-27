import sys
try:
    import psyco
    psyco.full()
except:
    print "psyco not running"
from cistematic.core import complement
from cistematic.genomes import Genome

delimiter = '|'
#delimiter = ':'

print "%s: version 3.5" % sys.argv[0]
if len(sys.argv) < 5:
    print "usage: python %s genome ucscModels outfilename maxBorder [-verbose] [-spacer num]\n" % sys.argv[0]
    print "\twhere spacer is by default 2, and maxBorder should be readlen - (2 * spacer)"
    print "\tdelimiter is set to %s - edit the code to change it, if necessary\n" % delimiter
    sys.exit(1)

genome = sys.argv[1]
datafilename = sys.argv[2]
outfilename = sys.argv[3]
# maxBorder should be readlen - 4
maxBorder = int(sys.argv[4])

doVerbose = False
if '-verbose' in sys.argv:
    doVerbose = True

spacer = 2
if '-spacer' in sys.argv:
    spacer = int(sys.argv[sys.argv.index('-spacer') + 1])
spacerseq = 'N' * spacer

datafile = open(datafilename)
#seqfile = open('knownGeneMrna.txt')
hg = Genome(genome)

spliceCountDict = {}
exonStartDict = {}
exonStopDict = {}
exonLengthDict = {}
nameToChromDict = {}
nameToComplementDict = {}
alreadySeen = {}
counter = 0

for line in datafile:
    fields = line.split()
    name = fields[0]
    spliceCount = int(fields[7]) - 1
    if spliceCount < 1:
        continue
    counter += spliceCount
    spliceCountDict[name] = spliceCount
    chrom = fields[1][3:]
    if chrom == 'chrM':
        continue
    nameToChromDict[name] = chrom
    if chrom not in alreadySeen:
        alreadySeen[chrom] = []
    nameToComplementDict[name] = fields[2]
    exonStarts = []
    exonStops = []
    for val in fields[8].split(',')[:-1]:
        exonStarts.append(int(val))
    for val in fields[9].split(',')[:-1]:
        exonStops.append(int(val))
    exonStartDict[name] = exonStarts
    exonStopDict[name] = exonStops
    exonLengths = []
    for index in range(spliceCount + 1):
        exonLengths.append(exonStops[index] - exonStarts[index])
    exonLengthDict[name] = exonLengths
print len(spliceCountDict)
print counter

missedCount = 0
depressedCount = 0
splicefileindex = 1
spliceCounter = 0
seenSpliceList = []
outfile = open(outfilename,'w')
for name in nameToChromDict:
    try:
        spliceCount = spliceCountDict[name]
    except:
        continue
    exonStarts = exonStartDict[name]
    exonStops = exonStopDict[name]
    exonLengths = exonLengthDict[name]
    chrom = nameToChromDict[name]
    exonCumulative = 0
    oldregionstart = 0
    for index in range(spliceCount):
        if (exonStops[index], exonStarts[index + 1]) in alreadySeen[chrom]:
            continue
        regionstart = exonStops[index] - maxBorder
        alreadySeen[chrom].append((exonStops[index], exonStarts[index + 1]))
        beforeLen = exonLengths[index]
        afterLen = exonLengths[index + 1]
        if (beforeLen + afterLen) < maxBorder + spacer:
            #print 'splice chr%s:%d-%d too short: %d' % (chrom, exonStops[index], exonStarts[index + 1], beforeLen + afterLen)
            missedCount += 1
            continue
        if (beforeLen + afterLen) < 2 * maxBorder:
            depressedCount += 1
        if beforeLen > maxBorder:
            beforeLen = maxBorder
        if afterLen > maxBorder:
            afterLen = maxBorder
        try:
            beforeSplice = hg.sequence(chrom, exonStops[index] - maxBorder, maxBorder)
            afterSplice = hg.sequence(chrom, exonStarts[index + 1], maxBorder)
            #beforeSplice = hg.sequence(chrom, exonStops[index] - beforeLen, beforeLen)
            #afterSplice = hg.sequence(chrom, exonStarts[index + 1], afterLen)
        except:
            if doVerbose:
                print "could not get chr%s:%d-%d" % (chrom, exonStops[index], exonStarts[index + 1])
            continue
        outstring = '>%s%s%d%s%d\n%s\n' % (name, delimiter, index, delimiter, regionstart, spacerseq + beforeSplice.upper() + afterSplice.upper() + spacerseq)
        outfile.write(outstring)
    splicefileindex += 1
    spliceCounter += 1
    if spliceCounter > 10000:
        print "%d genes" % splicefileindex
        spliceCounter = 0
    
outfile.close()

print "%d splices too short to be seen" % missedCount
print "%d splices will be under-reported" % depressedCount


