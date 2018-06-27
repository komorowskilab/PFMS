#
#  makewiggle.py
#  ENRAGE
#
import sys

print '%s: version 6.7' % sys.argv[0]

try:
    import psyco
    psyco.full()
except:
    print 'psyco not running'

from commoncode import *
from array import array

if len(sys.argv) < 4:
    print 'usage: python %s name rdsfile outfilename [-raw] [-color r,g,b] [-altcolor r,g,b] [-chrom chromID] [-shift bp] [-split] [-listfile filename] [-listprefix prefix] [-group groupName] [-startPriority prefix] [-skiprandom] [-nomulti] [-splices] [-singlebase] [-cache pages] [-enforceChr] [-stranded plus/minus/both] [-maxchunk Mbp]' % sys.argv[0]
    sys.exit(1)

name = sys.argv[1]
hitfilename = sys.argv[2]
outfilename = sys.argv[3]
colorString = ''
doNormalize = True
doSplit = False
doList = False
doSingle = False
withMulti = True
withSplices = False
listPrefix = ''
groupName=''
startPriority=0.01
priorityIncrement = 0.01
skipRandom = False
cachePages = -1
wigType = 'bedGraph'
#wigType = 'wiggle_0'
maxSpan = 20000000
shift = 0

if '-raw' in sys.argv:
    doNormalize = False

if '-split' in sys.argv:
    doSplit = True

if '-nomulti' in sys.argv:
    withMulti = False

if '-splices' in sys.argv:
    withSplices = True

if '-singlebase' in sys.argv:
    doSingle = True

if '-color' in sys.argv:
    colorField = sys.argv.index('-color') + 1
    try:
        colorString = ' color=' + sys.argv[colorField]
    except:
        pass

if '-altcolor' in sys.argv:
    colorField = sys.argv.index('-altcolor') + 1
    try:
        colorString += ' altcolor=' + sys.argv[colorField]
    except:
        pass

if '-listfile' in sys.argv:
    listfilename = sys.argv[sys.argv.index('-listfile') + 1]
    doList = True
if '-listprefix' in sys.argv:
    listPrefix = sys.argv[sys.argv.index('-listprefix') + 1]
    
chromLimit = False
if '-chrom' in sys.argv:
    chromField = sys.argv.index('-chrom') + 1
    chromLimit = True
    limitChrom = sys.argv[chromField]

if '-group' in sys.argv:
    groupName = 'group=' + sys.argv[sys.argv.index('-group') + 1]

if '-startPriority' in sys.argv:
    startPriority = float(sys.argv[sys.argv.index('-startPriority') + 1])

if '-skiprandom' in sys.argv:
    skipRandom = True

doCache = False
if '-cache' in sys.argv:
    doCache = True
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])

enforceChr = False
if '-enforceChr' in sys.argv:
    enforceChr = True

if '-maxchunk' in sys.argv:
    maxSpan = int(sys.argv[sys.argv.index('-maxchunk') + 1]) * 1000000

isStranded = False
strandedDirection = 'both'

if '-stranded' in sys.argv:
    isStranded = True
    nextArg =  sys.argv[sys.argv.index('-stranded') + 1]
    if nextArg == 'plus':
        strandedDirection = 'plusOnly'
    elif nextArg == 'minus':
        strandedDirection = 'minusOnly'
    print 'will keep track of %s strand(s)' % strandedDirection

if '-shift' in sys.argv:
    shift = int(sys.argv[sys.argv.index('-shift') + 1])
    print 'Will shift reads by +/- %d bp according to their sense' % shift
    name += 'shift=%d' % shift
    
hitRDS = readDataset(hitfilename, verbose = True, cache=doCache)

#sqlite default_cache_size is 2000 pages
if cachePages > hitRDS.getDefaultCacheSize():
    hitRDS.setDBcache(cachePages)

readlen = hitRDS.getReadSize()

if doNormalize:
    normalizeBy = len(hitRDS) / 1000000.
else:
    normalizeBy = 1.

if doList:
    listfile = open(listfilename, 'w')

priority = startPriority    
if not doSplit:
    outfile = open(outfilename,'w')
    if doList:
        listfile.write(listPrefix + outfilename + '\n')
    outfile.write('track type=%s name="%s" %s priority=%.3f visibility=full%s\n' % (wigType, name, groupName, priority, colorString)) 

chromList = hitRDS.getChromosomes()
chromList.sort()
for achrom in chromList:
    if enforceChr and ('chr' not in achrom):
        continue
    if chromLimit and achrom != limitChrom:
        continue
    if skipRandom and 'random' in achrom:
        continue
    if doSplit:
        outfile = open(outfilename + '.' + chrom,'w')
        if doList:
            listfile.write(listPrefix + outfilename + '.' + chrom + '\n')
        outfile.write('track type=%s name="%s %s" %s priority=%.3f visibility=full%s\n' % (wigType, name, chrom, groupName, priority, colorString))   
        priority += priorityIncrement  
    
    lastNT = hitRDS.getMaxCoordinate(achrom, doMulti=withMulti, doSplices=withSplices) + readlen
    spanStart = 0
    
    previousVal = 0
    previousStart = 1
    lineIndex = 0
    for spanStop in xrange(maxSpan, lastNT+maxSpan, maxSpan):
        if spanStop > lastNT:
            spanStop = lastNT
        print achrom, spanStart, spanStop
        chromModel = hitRDS.getChromProfile(achrom, spanStart, spanStop, withMulti, withSplices, normalizeBy, isStranded, strandedDirection, shiftValue=shift)
    
        for index in xrange(len(chromModel)):
            currentVal = chromModel[index]
            if doSingle:
                outline = '%s %d %.4f\n' % (achrom, spanStart + index, currentVal)
                outfile.write(outline)
                continue
            if currentVal == previousVal:
                continue
            if currentVal != previousVal:
                if previousVal != 0:
                    lastpos = index + spanStart
                    outline = '%s %d %d %.4f\n' % (achrom, previousStart, lastpos, previousVal)
                    outfile.write(outline)
                    lineIndex += 1
                previousVal = currentVal
                previousStart = index + spanStart
        currentVal = 0
        del chromModel
        spanStart = spanStop + 1
    if doSplit:
        outfile.close()
    if doSingle:
        print index + 1
    else:
        print lineIndex
if not doSplit:
    outfile.close()
if doList:
    listfile.close()

