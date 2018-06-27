import sys
try:
    import psyco
    psyco.full()
except:
	print 'psyco not running'
from cistematic.core.motif import Motif, hasMotifExtension
from cistematic.core import complement
from cistematic.genomes import Genome
from commoncode import *

print '%s: version 2.4' % sys.argv[0]
if len(sys.argv) < 6:
    print 'usage: python %s genome motifFile motThreshold regionfile siteOutfile [-dataset chipRDS] [-cache] [-best] [-usepeak] [-rank] [-printseq] [-nomerge] [-markov1] ' % sys.argv[0]
    sys.exit(1)

genome = sys.argv[1]
motfilename = sys.argv[2]
motThreshold = float(sys.argv[3])

doMarkov1 = False
if '-markov1' in sys.argv:
    doMarkov1 = True

if motThreshold < 1.0 and doMarkov1:
    print 'motThreshold should be between 1.0 and 10.0 for markov1'
    sys.exit(1)
elif motThreshold < 55.0 and not doMarkov1:
    print 'motThreshold should be between 55 and 99 for a regular PSFM'
    sys.exit(1)

#pvalfilename = 'pvalueunamp.txt'
infilename = sys.argv[4]
outfilename = sys.argv[5]

bestOnly = False
if '-best' in sys.argv:
    print "will only report the best position for each region"
    bestOnly = True

usePeak = False
if '-usepeak' in sys.argv:
    print "will use peak position and height from regions file"
    usePeak = True

useRank = False
if '-rank' in sys.argv:
    if usePeak:
        print "will return region ranking based on peak height ranking"
        useRank = True
    else:
        print "ignoring '-rank': can only use ranking when using a region file with peak position and height"

doRDS = False
if '-dataset' in sys.argv:
    chipfilename = sys.argv[sys.argv.index('-dataset') + 1]
    doRDS = True

doCache = False
if '-cache' in sys.argv:
    doCache = True

printSeq = False
if '-printseq' in sys.argv:
    printSeq = True
maxPvalue = 0.0001

mot = Motif('',motifFile = motfilename)
motLen = len(mot)
bestScore = mot.bestConsensusScore()

if hasMotifExtension:
    print "will use cistematic.core.motif C-extension to speed up motif search"

hg = Genome(genome)

# minHits=-1 will force regions to be used regardless
# maxDist= 0 prevents merging of non-overlapping regions
if '-nomerge' in sys.argv:
    regions = getMergedRegions(infilename, maxDist=0, minHits=-1, verbose=True, doMerge=False, keepPeak=usePeak)
else:
    regions = getMergedRegions(infilename, maxDist=0, minHits=-1, verbose=True, keepPeak=usePeak)

if doRDS:
    hitRDS = readDataset(chipfilename, verbose = True, cache=doCache)

outfile = open(outfilename,'w')

index = 0
regionList = []

for chrom in regions:
    if 'rand' in chrom or 'M' in chrom:
        continue
    if usePeak:
        for (start, stop, length, peakPos, peakHeight) in regions[chrom]:
            regionList.append((peakHeight, chrom, start, length, peakPos))
    else:
        for (start, stop, length) in regions[chrom]:
            regionList.append((chrom, start, length))

if usePeak:
    regionList.sort()
    regionList.reverse()
notFoundIndex = 0
currentChrom = ''
count = 0
for tuple in regionList:
    if usePeak:
        (rpeakheight, rchrom, start, length, rpeakpos) = tuple
    else:
        (rchrom, start, length) = tuple
    try:
        seq = hg.sequence(rchrom, start, length)
    except:
        print "couldn't retrieve %s %d %d - skipping" % (rchrom, start, length)
        continue
    count += 1
    numHits = -1
    if usePeak:
        peakpos = rpeakpos
        if useRank:
            numHits = count
        else:
            numHits = rpeakheight
    elif doRDS:
        if rchrom != currentChrom:
            fullchrom = 'chr' + rchrom
            hitDict = hitRDS.getReadsDict(chrom=fullchrom)
            currentChrom = rchrom
        (topPos, numHits, smoothArray, numPlus) = findPeak(hitDict[rchrom], start, length)
        if len(topPos) == 0:
            print 'topPos error'
        peakpos = topPos[0]
    found = []
    if doMarkov1:
        matches = mot.locateMarkov1(seq, motThreshold)
    else:
        matches = mot.locateMotif(seq, motThreshold)
    for (pos, sense) in matches:
        alreadyFound = False
        for (fpos, fdist) in found:
            if pos + start == fpos:
                alreadyFound = True
        if not alreadyFound:
            if usePeak:
                found.append((start + pos, start + pos  + motLen/2 - peakpos))
            elif doRDS:
                found.append((start + pos, pos  + motLen/2 - peakpos))
            else:
                found.append((start + pos, -1))
                            
    foundValue = False
    bestList = []
    for (foundpos, peakdist) in found:
        seq = hg.sequence(rchrom, foundpos, motLen)
        foundValue = True
        (front, back) = mot.scoreMotif(seq)
        sense = '+'
        if front >= back:
            score = int(100 * front / bestScore)
        else:
            score = int(100 * back / bestScore)
            sense = '-'
            #foundpos += 1
            seq = complement(seq)
        if printSeq:
            print seq
        outline = 'chr%s:%d-%d\t%d\t%d\t%d\tchr%s:%d-%d\t%s\n' % (rchrom, foundpos, foundpos + motLen - 1, score, numHits, peakdist, rchrom, start, start + length, sense)
        if bestOnly:
            bestList.append((abs(peakdist), outline))
        else:
            outfile.write(outline)
    if bestOnly and foundValue:
        bestList.sort()
        outfile.write(bestList[0][1])
    if not foundValue:
        if printSeq:
            print 'could not find a %s site for %s:%d-%d' % (mot.tagID, rchrom, start, start+ length)
        notFoundIndex += 1
    if (count % 10000) == 0 and not printSeq:
        print count
outfile.close()
print 'did not find motif in %d regions' % notFoundIndex
