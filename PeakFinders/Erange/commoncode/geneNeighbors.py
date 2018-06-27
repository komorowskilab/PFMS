#
#  geneNeighbors.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    pass

import sys
print '%s: version 2.4' % sys.argv[0]
from commoncode import *

if len(sys.argv) < 3:
    print 'usage: python %s genome outfilename [-regions acceptfile] [-downstream bp] [-upstream bp] [-mindist bp] [-minlocus bp] [-maxlocus bp] [-samesense]' % sys.argv[0]
    sys.exit(1)

from cistematic.genomes import Genome
from cistematic.core.geneinfo import geneinfoDB

genome = sys.argv[1]
outfilename = sys.argv[2]

acceptDict = {}
if '-regions' in sys.argv:
    acceptField = sys.argv.index('-regions') + 1
    acceptfile = sys.argv[acceptField]
    acceptDict = getMergedRegions(acceptfile, maxDist=0, keepLabel=True, verbose=True)

checkSense = False
if '-samesense' in sys.argv:
    checkSense = True

downMax = 10000000
if '-downstream' in sys.argv:
    downField = sys.argv.index('-downstream') + 1
    downMax = int(sys.argv[downField])

upMax = 10000000
if '-upstream' in sys.argv:
    upField = sys.argv.index('-upstream') + 1
    upMax = int(sys.argv[upField])

minDist = 0
if '-mindist' in sys.argv:
    minField = sys.argv.index('-mindist') + 1
    minDist = int(sys.argv[minField])

minLocus = -1
if '-minlocus' in sys.argv:
    minlfield = sys.argv.index('-minlocus') + 1
    minLocus = int(sys.argv[minlfield])

maxLocus = 10000000
if '-maxlocus' in sys.argv:
    maxlfield = sys.argv.index('-maxlocus') + 1
    maxLocus = int(sys.argv[maxlfield])

hg = Genome(genome)
idb = geneinfoDB(cache=True)

geneinfoDict = idb.getallGeneInfo(genome)
locusByChromDict = getLocusByChromDict(hg, additionalRegionsDict=acceptDict, keepSense=True)

gidList = hg.allGIDs()
gidList.sort()
for chrom in acceptDict:
    for (label, start, stop, length) in acceptDict[chrom]:
        if label not in gidList:
            gidList.append(label)

index = 0
outfile = open(outfilename,'w')
chromList = locusByChromDict.keys()
chromList.sort()
for chrom in chromList:
    if len(locusByChromDict[chrom]) < 3 or 'NT' in chrom or 'MT' in chrom:
        continue
    print chrom + ' ',
    
    prevStop = locusByChromDict[chrom][0][1]
    prevGID = locusByChromDict[chrom][0][2]
    if 'FAR' not in prevGID:
        symbol = 'LOC' + prevGID
        geneinfo = ''
        try:
            geneinfo = geneinfoDict[prevGID]
            symbol = geneinfo[0][0]
        except:
            #print 'problem with %s' % (gid)
            pass
    else:
        symbol = prevGID
    prevGID = symbol
    prevSense = locusByChromDict[chrom][0][4]
    
    currentStart = locusByChromDict[chrom][1][0]
    currentStop = locusByChromDict[chrom][1][1]
    currentGID = locusByChromDict[chrom][1][2]
    if 'FAR' not in currentGID:
        symbol = 'LOC' + currentGID
        geneinfo = ''
        try:
            geneinfo = geneinfoDict[currentGID]
            symbol = geneinfo[0][0]
        except:
            #print 'problem with %s' % (gid)
            pass
    else:
        symbol = currentGID
    currentGID = symbol
    currentGlen = locusByChromDict[chrom][1][3]
    currentSense = locusByChromDict[chrom][1][4] 
    
    for (nextStart, nextStop, nextGID, nextGlen, nextSense) in locusByChromDict[chrom][2:]:
        if 'FAR' not in nextGID:
            symbol = 'LOC' + nextGID
            geneinfo = ''
            try:
                geneinfo = geneinfoDict[nextGID]
                symbol = geneinfo[0][0]
            except:
                #print 'problem with %s' % (gid)
                pass
        else:
            symbol = nextGID
        nextGID = symbol
        leftDist = currentStart - prevStop
        rightDist = nextStart - currentStop
        if (currentSense == 'F' and minDist < leftDist < upMax and minDist < rightDist < downMax) or (currentSense == 'R' and minDist < rightDist < upMax and minDist < leftDist < downMax):
            if not checkSense or currentSense == nextSense:
                if minLocus <= currentGlen <= maxLocus:
                    outfile.write('%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\n' % (currentGID, currentSense, prevGID, prevSense, leftDist, nextGID, nextSense, rightDist))
                    index += 1
        prevStop = currentStop
        prevGID = currentGID
        prevSense = currentSense
        currentStart = nextStart
        currentStop = nextStop
        currentGID = nextGID
        currentGlen = nextGlen
        currentSense = nextSense
outfile.close()

print '\n%d genes matched' % index