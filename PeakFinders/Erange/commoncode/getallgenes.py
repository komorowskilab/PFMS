try:
	import psyco
	psyco.full()
except:
	print 'psyco not running'
from cistematic.core import genesIntersecting, featuresIntersecting, cacheGeneDB, uncacheGeneDB
from cistematic.core.geneinfo import geneinfoDB
from cistematic.genomes import Genome
import sys

print '%s: version 5.5' % sys.argv[0]
if len(sys.argv) < 4:
    print 'usage: python %s genome regionfile outfile [-radius bp] [-nomatch nomatchfile] -trackfar -stranded -cache -compact [-step dist] [-startField colID]' % sys.argv[0]
    sys.exit(1)

nomatchfilename = ''
genome = sys.argv[1]
infilename = sys.argv[2]
outfilename = sys.argv[3]
maxRadius = 20002
doCache = False
trackFar = False
trackStrand = False

if '-trackfar' in sys.argv:
    trackFar = True
    print 'trackfar = True'

if '-stranded' in sys.argv:
    trackStrand = True
    
compact=False
if '-compact' in sys.argv:
    compact = True

colID = 1
if '-startField' in sys.argv:
    colID = int(sys.argv[sys.argv.index('-startField') + 1])
    
if '-radius' in sys.argv:
    try:
        maxRadius = int(sys.argv[sys.argv.index('-radius') + 1])
        print "using a radius of %d" % maxRadius
    except:
        pass

step = maxRadius - 2
if '-step' in sys.argv:
    step = int(sys.argv[sys.argv.index('-step') + 1])

if '-nomatch' in sys.argv:
        nomatchfilename = sys.argv[sys.argv.index('-nomatch') + 1]

extendGenome = ''
replaceModels = False
if '-models' in sys.argv:
    extendGenome = sys.argv[sys.argv.index('-models') + 1]
    if '-replacemodels' in sys.argv:
        replaceModels = True
        print "will replace gene models with %s" % extendGenome
    else:
        print "will extend gene models with %s" % extendGenome

if '-cache' in sys.argv:
    doCache = True
    idb = geneinfoDB(cache=True)
    #print "cached %s" % genome
else:
    idb = geneinfoDB()

infile = open(infilename)
outfile = open(outfilename,'w')

if genome == 'dmelanogaster':
    geneinfoDict = idb.getallGeneInfo(genome, infoKey='locus')
else:
    geneinfoDict = idb.getallGeneInfo(genome)

posList = []
altPosDict = {}
altPosRevDict = {}
posLine = {}
posStrand = {}
altPosList = []

index = 0
for line in infile:
    if line[0] == '#':
        continue
    fields = line.split('\t')
    if compact:
        (chrom, pos) = fields[colID].split(':')
        chrom = chrom[3:]
        (start, stop) = pos.split('-')
        pos = (chrom, int(start))
        altPos = (chrom, int(stop))
    else:
        try:
            chrom = fields[colID][3:]
        except:
            print line
            continue
        pos = (chrom, int(fields[colID + 1]))
        altPos = (chrom, int(fields[colID + 2]))
    
    altPosDict[pos] = altPos
    altPosRevDict[altPos] = pos
    posList.append(pos)
    posList.append(altPos)
    altPosList.append(altPos)
    posLine[pos] = line
    if trackStrand:
        if 'RNAFARP' in line:
            posStrand[pos] = '+'
            posStrand[altPos] = '+'
        else:
            posStrand[pos] = '-'
            posStrand[altPos] = '-'
	
geneList = []
geneDict = {}
if maxRadius < step:
    step = maxRadius - 2

hg = Genome(genome, inRAM=True)
if extendGenome != '':
    hg.extendFeatures(extendGenome, replace = replaceModels)
    
geneannotDict = hg.allAnnotInfo()
#featureTypes = ['CDS'] + hg.getFeatureTypes('UT%')
featureTypes = ['CDS', 'UTR']
for radius in range(1, maxRadius, step):
    print 'radius %d' % radius
    print len(posList)
    if radius == 1:
            posDict = genesIntersecting(genome, posList, extendGen=extendGenome, replaceMod=replaceModels)
    else:
            posDict = featuresIntersecting(genome, posList, radius, 'CDS', extendGen=extendGenome, replaceMod=replaceModels) 
            posDict2 = featuresIntersecting(genome, posList, radius, 'UTR', extendGen=extendGenome, replaceMod=replaceModels)
            for apos in posDict2:
                try: 
                    posDict[apos] += posDict2[apos]
                    posDict[apos].sort()
                except:
                    posDict[apos] = posDict2[apos]
    for pos in posDict:
        geneID  = ''
        if len(posDict[pos]) == 1:
            if trackStrand:
                if posStrand[pos] == posDict[pos][0][-1]:
                    geneID = posDict[pos][0][0]
            else:
                geneID = posDict[pos][0][0]
        elif len(posDict[pos]) > 1 and not trackStrand:
            (chrom, loc) = pos
            bestres = posDict[pos][0]
            dist1 = abs(bestres[3] - loc)
            dist2 = abs(bestres[4] - loc)
            if dist1 < dist2:
                bestdist = dist1
            else:
                bestdist = dist2
            for testres in posDict[pos]:
                testdist1 = abs(testres[3] - loc)
                testdist2 = abs(testres[4] - loc)
                if testdist1 < testdist2:
                    testdist = testdist1
                else:
                    testdist = testdist2
                if testdist < bestdist:
                    bestdist = testdist
                    bestres = testres
            geneID = bestres[0]
        elif len(posDict[pos]) > 1:
            (chrom, loc) = pos
            bestres = posDict[pos][0]
            dist1 = abs(bestres[3] - loc)
            dist2 = abs(bestres[4] - loc)
            bestStrand = posDict[pos][-1]
            if dist1 < dist2:
                bestdist = dist1
            else:
                bestdist = dist2
            for testres in posDict[pos]:
                testdist1 = abs(testres[3] - loc)
                testdist2 = abs(testres[4] - loc)
                testStrand = testres[-1]
                if testdist1 < testdist2:
                    testdist = testdist1
                else:
                    testdist = testdist2
                if bestStrand != posStrand[pos] and testStrand == posStrand[pos]:
                    bestdist = testdist
                    bestres = testres
                    bestStrand = testStrand
                elif testdist < bestdist:
                    bestdist = testdist
                    bestres = testres
            if bestStrand == posStrand[pos]:
                geneID = bestres[0]
        
        if geneID != '':
            try:
                if genome == 'dmelanogaster':
                    symbol = geneinfoDict['Dmel_' + geneID][0][0]
                else:
                    symbol = geneinfoDict[geneID][0][0]
            except:
                try:
                    symbol = geneannotDict[(genome, geneID)][0]
                except:
                    symbol = 'LOC' + geneID
        else:
            continue
        
        if pos in altPosList and pos in posList:
            posList.remove(pos)
            if pos not in altPosRevDict:
                continue
            if altPosRevDict[pos] in posList:
                posList.remove(altPosRevDict[pos])
            pos = altPosRevDict[pos]
        elif pos in posList:
            posList.remove(pos)
            if pos not in altPosDict:
                print pos
                continue
            if altPosDict[pos] in posList:
                posList.remove(altPosDict[pos])
        else:
            continue
            
        if (symbol, geneID) not in geneList:
            geneList.append((symbol, geneID))
            geneDict[(symbol, geneID)] = []
        if pos not in geneDict[(symbol, geneID)]:
            geneDict[(symbol, geneID)].append(pos)

for (symbol, geneID) in geneList:
    geneDict[(symbol, geneID)].sort()
    seenLine = []
    for pos in geneDict[(symbol, geneID)]:
        if pos in altPosRevDict:
            pos = altPosRevDict[pos]
        if posLine[pos] in seenLine:
            continue
        if '\t' in symbol:
            symbol = symbol.replace('\t','|')
        if ' ' in symbol:
            symbol = symbol.replace(' ','_')
        line = '%s %s %s' % (symbol, geneID, posLine[pos])
        seenLine.append(posLine[pos])
        outfile.write(line)
    #print line

matchIndex = 0
if nomatchfilename != '':
    nomatchfile = open(nomatchfilename, 'w')
prevStart = 0
prevChrom = ''
farIndex = 0
start = 0
for pos in posList:
    if pos not in altPosList:
        if nomatchfilename != '':
            nomatchfile.write(posLine[pos])
        matchIndex += 1
        # need to add strand tracking here.....
        if trackFar:
            (chrom, start) = pos
            if chrom != prevChrom:
                farIndex += 1
                prevChrom = chrom
            elif abs(int(start) - prevStart) > maxRadius:
                farIndex += 1
            line = 'FAR%d %d %s' % (farIndex, -1 * farIndex, posLine[pos])
            outfile.write(line)
        prevStart = int(start)
if nomatchfilename != '':
    nomatchfile.close()
print '%d sites without a gene within radius of %d' % (matchIndex, radius)

