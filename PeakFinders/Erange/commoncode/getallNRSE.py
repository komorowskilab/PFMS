import sys

try:
    import psyco
    psyco.full()
except:
    print 'psyco not running'
from cistematic.core import complement
from cistematic.core.motif import Motif
from cistematic.genomes import Genome
from commoncode import *

print '%s: version 3.4' % sys.argv[0]
if len(sys.argv) < 4:
    print 'usage: python %s genome regionfile siteOutfile [-dataset chipRDS] [-min heightRPM] [-minfraction heightfraction] [-plot badregionPlot.png] [-cache] [-raw] [-verbose] [-markov1] [-peakdist bp] [-fullOnly] [-motifdir directory]' % sys.argv[0]
    sys.exit(1)

doPlot = False
if '-plot' in sys.argv:
    import matplotlib
    matplotlib.use('Agg')
    
    from pylab import *
    doPlot = True
    plotname = sys.argv[sys.argv.index('-plot') + 1]

genome = sys.argv[1]
#pvalfilename = 'pvalueunamp.txt'
infilename = sys.argv[2]
outfilename = sys.argv[3]

doCache = False
if '-cache' in sys.argv:
    doCache = True
    
normalize = True
normalizeBy = 1.
if '-raw' in sys.argv:
    normalize = False

doVerbose = False
if '-verbose' in sys.argv:
    doVerbose = True
    
motifDir = './'
if '-motifdir' in sys.argv:
    motifDir = sys.argv[sys.argv.index('-motifdir') + 1]
    if motifDir[-1] != '/':
        motifDir += '/'

doDataset = False
chipfilename = ''
if '-dataset' in sys.argv:
    chipfilename = sys.argv[sys.argv.index('-dataset') + 1]
    hitRDS = readDataset(chipfilename, verbose = doVerbose, cache=doCache)
    doDataset = True
    if normalize:
        normalizeBy = len(hitRDS) / 1000000.

minHeight = -2.
if '-min' in sys.argv:
    minHeight = float(sys.argv[sys.argv.index('-min') + 1])

minFraction = -2.
if '-minfraction' in sys.argv:
    minFraction = float(sys.argv[sys.argv.index('-minfraction') + 1])
    if minFraction > 1.:
        minFraction /= 100.
        print 'scaling minFraction to %.2f' % minFraction

doMarkov1 = False
if '-markov1' in sys.argv:
    doMarkov1 = True

maxpeakdist = 101
enforcePeakDist = False
if '-peakdist' in sys.argv:
    enforcePeakDist = True
    maxpeakdist = int(sys.argv[sys.argv.index('-peakdist') + 1])

fullOnly = False
if '-fullOnly' in sys.argv:
    fullOnly = True
    
#mot = Motif('',motifFile = motifDir + 'NRSE2.mot')
#motL = Motif('',motifFile = motifDir + 'NRSE2left.mot')
#motR = Motif('',motifFile = motifDir + 'NRSE2right.mot')
mot = Motif('',motifFile = motifDir + 'NRSE3.mot')
motL = Motif('',motifFile = motifDir + 'NRSE3left.mot')
motR = Motif('',motifFile = motifDir + 'NRSE3right.mot')
bestScore = mot.bestConsensusScore()
bestLeft = motL.bestConsensusScore()
bestRight = motR.bestConsensusScore()

hg = Genome(genome)

regions = getMergedRegions(infilename, maxDist=0, minHits=-1, verbose=doVerbose, doMerge=False)

outfile = open(outfilename,'w')
outfile.write('#dataset: %s\tregions:%s\tnormalize: %s\tmarkov1: %s\n' % (chipfilename, infilename, normalize, doMarkov1))
outfile.write('#enforcePeakDist: %s\tpeakdist: %d bp\tfullOnly: %d bp\n' % (enforcePeakDist, maxpeakdist, fullOnly))
outfile.write('#site\tscore\tleftscore\trightscore\tRPM\tpeakDist\ttype\theight\tfractionHeight\tregion\tsense\tseq\n')
		
countList = []
posList = []

index = 0
regionList = []

for rchrom in regions:
    if 'rand' in rchrom or 'M' in rchrom or 'hap' in rchrom:
        continue
    for (start, stop, length) in regions[rchrom]:
        regionList.append((rchrom, start, length))
	
notFoundIndex = 0
currentChrom = ''
for (rchrom, start, length) in regionList:
    seq = hg.sequence(rchrom, start, length)
    if doDataset:
        if rchrom != currentChrom:
            fullchrom = 'chr' + rchrom
            hitDict = hitRDS.getReadsDict(chrom=fullchrom, withWeight=True, doMulti=True)
            currentChrom = rchrom
        (topPos, numHits, smoothArray, numPlus) = findPeak(hitDict[rchrom], start, length, doWeight=True)
        if len(topPos) == 0:
            print 'topPos error'
        peakpos = topPos[0]
        peakscore = smoothArray[peakpos]
        if peakscore == 0.:
            peakscore = -1.
        if normalize:
            numHits /= normalizeBy
            peakscore /= normalizeBy
    else:
        peakpos = length
        peakscore = -1
        numHits = 0
        smoothArray = [0.] * length
    found = []
    if doMarkov1:
        lefts = motL.locateMarkov1(seq, 3.)
        rights = motR.locateMarkov1(seq, 3.)
    else:
        lefts = motL.locateMotif(seq, 70)
        rights = motR.locateMotif(seq, 70)
    allhalfs = [(v0, v1, 'L') for (v0, v1) in lefts] + [(v0, v1, 'R') for (v0, v1) in rights]
    allhalfs.sort()
	
    # look for canonicals and non-canonicals
    if len(allhalfs) > 1:
        (firstpos, firstsense, firsttype) = allhalfs[0]
        for (secondpos, secondsense, secondtype) in allhalfs[1:]:
            if enforcePeakDist:
                withinDistance = False
                for aPos in topPos:
                    if abs(firstpos - aPos) < maxpeakdist or abs(secondpos - aPos) < maxpeakdist:
                        withinDistance = True
                if not withinDistance:
                    firstpos = secondpos
                    firstsense = secondsense
                    firsttype = secondtype
                    continue
            if firsttype == 'L':
                dist = secondpos - firstpos + 2
            else:
                dist = secondpos - firstpos -1
            if firstsense == secondsense and dist in [9, 10, 11, 16, 17, 18, 19]:
                if (firsttype == 'L' and secondtype == 'R' and secondsense == 'F'):
                    found.append((start + firstpos, firstpos - peakpos + (dist + 10)/2, dist))
                if (firsttype == 'R' and secondtype == 'L' and secondsense == 'R'):
                    found.append((start + firstpos, firstpos  - peakpos + (dist + 10)/2, dist))
            firstpos = secondpos
            firstsense = secondsense
            firsttype = secondtype
		
    # did we miss any 70%+ matches ?
    if doMarkov1:
        matches = mot.locateMarkov1(seq, 3.5)
    else:
        matches = mot.locateMotif(seq, 70)
    for (pos, sense) in matches:
        alreadyFound = False
        for (fpos, fpeakdist, fdist) in found:
            if pos + start == fpos:
                alreadyFound = True
        if not alreadyFound:
            if enforcePeakDist:
                withinDistance = False
                for aPos in topPos:
                    if abs(firstpos - aPos) < maxpeakdist or abs(secondpos - aPos) < maxpeakdist:
                        withinDistance = True
                        thePos = aPos
                if withinDistance:
                    found.append((start + pos, pos - thePos + 10, 11))

            else:
                found.append((start + pos, pos - peakpos + 10, 11))
	
    # we'll now accept half-sites within maxpeakdist bp of peak if using a dataset, else all
    if len(found) == 0 and not fullOnly:
        bestone = -1
        if not doDataset:
            bestdist = maxpeakdist
        else:
            bestdist = length
        index = 0
        for (pos, sense, type) in allhalfs:
            if doDataset:
                for aPos in topPos:
                    if abs(pos - aPos) < bestdist:
                        bestdist = abs(pos - aPos)
                        bestone = index
                        peakpos = aPos
            else:
                found.append((start + allhalfs[index][0], allhalfs[index][0] + 5 - peakpos, 0))
                    
            index += 1
        if (doDataset and bestdist < 101):
            try:
                found.append((start + allhalfs[bestone][0], allhalfs[bestone][0] + 5 - peakpos, 0))
            except:
                continue
    # see if we found an acceptable match
    foundValue = False
    for (foundpos, posdist, dist) in found:
        # get a score for 21-mer, report
        seq = hg.sequence(rchrom, foundpos, 21)
        # height will be measured from the center of the motif
        height = -2.
        for pos in range(10 + dist):
            try:
                currentHeight = smoothArray[int(peakpos + posdist + pos)]
            except: 
                pass
            if currentHeight > height:
                height = currentHeight
        if normalize:
            height /= normalizeBy
        fractionHeight = height / peakscore
        if height < minHeight or fractionHeight < minFraction:
            continue
        foundValue = True
        (front, back) = mot.scoreMotif(seq)
        sense = '+'
        if front > back:
            score = int(100 * front / bestScore)
            theseq = hg.sequence(rchrom, foundpos, 10 + dist)
        else:
            score = int(100 * back / bestScore)
            theseq = complement(hg.sequence(rchrom, foundpos, 10 + dist))
            sense = '-'
            foundpos + 1
        leftScore = -1.
        rightScore = -1.
        leftseq = ''
        rightseq = ''
        if dist > 0:
            testseq = hg.sequence(rchrom, foundpos, 10 + dist)
            if sense == '-':
                testseq = complement(testseq)
            leftseq = testseq[:9]
            rightseq = testseq[dist-2:]
        elif dist == 0:
            testseq = hg.sequence(rchrom, foundpos, 12)
            if sense == '-':
                testseq = complement(testseq)
                leftseq = testseq[3:]
            else:
                leftseq = testseq[:9]
            rightseq = testseq
        (lfront, lback) = motL.scoreMotif(leftseq)
        #print (leftseq, lfront, lback)
        (rfront, rback) = motR.scoreMotif(rightseq)
        if lfront > lback:
            leftScore = int(100 * lfront) / bestLeft
            leftSense = '+'
        else:
            leftScore = int(100 * lback) / bestLeft
            leftSense = '-'
        if rfront > rback:
            rightScore = int(100 * rfront) / bestRight
            rightSense = '+'
        else:
            rightScore = int(100 * rback) / bestRight
            rightSense = '-'
        if dist != 11:
            if rightScore > leftScore:
                sense = rightSense
            else:
                sense = leftSense
            if sense == '-' and dist > 0:
                theseq = complement(hg.sequence(rchrom, foundpos, 10 + dist))
        
        outline = 'chr%s:%d-%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\tchr%s:%d-%d\t%s\t%s' % (rchrom, foundpos, foundpos + 9 + dist, score, leftScore, rightScore, numHits, posdist, dist, height, fractionHeight, rchrom, start, start + length, sense, theseq)
        if doVerbose:
            print outline
        outfile.write(outline + '\n')
	
    # we didn't find a site - draw region
    if not foundValue and doVerbose:
        outline = '#no predictions for %s:%d-%d %d %.2f' % (rchrom, start, start + length, numHits, peakscore)
        print outline
        outfile.write(outline + '\n')
    if not foundValue and doPlot:
        drawarray = [val + notFoundIndex for val in smoothArray]
        drawpos = [drawarray[val] for val in topPos]
        plot(drawarray,'b')
        plot(topPos, drawpos,'r.')
        goodmatches = mot.locateMotif(seq, 75)
        if len(goodmatches) > 0:
            print topPos
            print goodmatches
            drawgood = []
            drawgoody = []
            for (mstart, sense) in goodmatches:
                drawgood.append(mstart)
                drawgoody.append(drawarray[mstart])
            plot(drawgood, drawgoody, 'g.')
        notFoundIndex -= 30

outfile.close()
if doPlot:
    savefig(plotname)
