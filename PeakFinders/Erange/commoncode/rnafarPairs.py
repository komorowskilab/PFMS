#
#  RNAFARpairs.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 11/2/08.
#

try:
    import psyco
    psyco.full()
except:
    pass

from commoncode import *
from cistematic.core.geneinfo import geneinfoDB
from cistematic.genomes import Genome
import sys, time

print '%s: version 3.6' % sys.argv[0]

if len(sys.argv) < 2:
    print 'usage: python %s genome goodfile rdsfile outfile [-verbose] [-cache]' % sys.argv[0]
    print '\tlooks at all chromosomes simultaneously: is both slow and takes up large amount of RAM'
    sys.exit(1)

minDist = 0
genome = sys.argv[1]
goodfilename = sys.argv[2]
rdsfile = sys.argv[3]
outfilename = sys.argv[4]

sameChromOnly = False
if '-sameChrom' in sys.argv:
    sameChromOnly = True

doVerbose = False
if '-verbose' in sys.argv:
    doVerbose = True

doCache = False
if '-cache' in sys.argv:
    doCache = True

maxDist = 500000
if '-maxDist' in sys.argv:
    maxDist = int(sys.argv[sys.argv.index('-maxDist') + 1])

goodDict = {}
goodfile = open(goodfilename)
for line in goodfile:
    fields = line.split()
    goodDict[fields[0]] = line

RDS = readDataset(rdsfile, verbose = True, cache=doCache)
rdsChromList = RDS.getChromosomes()

if doVerbose:
    print time.ctime()

distinct = 0
total = 0
outfile = open(outfilename,'w')

idb = geneinfoDB()
if genome == 'dmelanogaster':
    geneinfoDict = idb.getallGeneInfo(genome, infoKey='locus')
else:
    geneinfoDict = idb.getallGeneInfo(genome)
hg = Genome(genome)
geneannotDict = hg.allAnnotInfo()

assigned = {}
farConnected = {}
for achrom in rdsChromList:
    if achrom == 'chrM':
        continue
    print achrom
    uniqDict = RDS.getReadsDict(fullChrom=True, chrom=achrom, noSense=True, withFlag=True, withPairID=True, doUniqs=True, readIDDict=True)
    if doVerbose:
        print len(uniqDict), time.ctime()    
    
    for readID in uniqDict:
        readList = uniqDict[readID]
        if len(readList) == 2:
            total += 1
            (start1, flag1, pair1) = readList[0]
            (start2, flag2, pair2) = readList[1]
            
            if flag1 != flag2:
                dist = abs(start1 - start2)
                if flag1 != 'NM' and flag2 != 'NM' and dist < maxDist:
                    geneID = ''
                    saw1 = False
                    saw2 = False
                    if flag1 in goodDict:
                        geneID = flag2
                        farFlag = flag1
                        saw1 = True
                    if flag2 in goodDict:
                        geneID = flag1
                        farFlag = flag2
                        saw2 = True
                    if saw1 or saw2:
                        total += 1
                    if saw1 and saw2:
                        if flag1 < flag2:
                            geneID = flag1
                            farFlag = flag2
                        else:
                            geneID = flag2
                            farFlag = flag1
                        if geneID in farConnected:
                            farConnected[geneID].append(farFlag)
                        else:
                            farConnected[geneID] = [farFlag]
                    elif geneID != '':
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
                        symbol = symbol.strip()
                        symbol = symbol.replace(' ','|')
                        symbol = symbol.replace('\t','|')
                        if farFlag not in assigned:
                            assigned[farFlag] = (symbol, geneID)
                            print '%s %s %s' % (symbol, geneID, goodDict[farFlag].strip())
                            outfile.write('%s %s %s' % (symbol, geneID, goodDict[farFlag]))
                            distinct += 1

farIndex = 0
for farFlag in farConnected:
    geneID = ''
    symbol = ''
    idList = [farFlag] + farConnected[farFlag]
    theID = False
    for oneID in idList:
        if oneID in assigned:
            (symbol, geneID) = assigned[oneID]
    if geneID == '':
        farIndex += 1
        symbol = 'FAR%d' % farIndex
        geneID = -1 * farIndex
    for oneID in idList:
        if oneID not in assigned:
            print '%s %s %s' % (symbol, geneID, goodDict[oneID].strip())
            outfile.write('%s %s %s' % (symbol, geneID, goodDict[oneID]))
            distinct += 1
            assigned[oneID] = (symbol, geneID)
        
                
for farFlag in goodDict:
    if farFlag not in assigned:
        farIndex += 1
        line = 'FAR%d %d %s' % (farIndex, -1 * farIndex, goodDict[farFlag])
        print line.strip()
        outfile.write(line)


outfile.write("#distinct: %d\ttotal: %d\n" % (distinct, total))
#total = float(total)
#if total < 1:
#    total = 1.
#outfile.write('#distal %2.2f pct\tdiffChrom %2.2f pct\n' % ((100. * distal/total), (100. * diffChrom/total)))
outfile.close()
print "distinct: %d\ttotal: %d" % (distinct, total)
#print "distal: %2.2f pct\tdiffChrom: %2.2f pct\n" % ((100. * distal/total), (100. * diffChrom/total))
print time.ctime()
