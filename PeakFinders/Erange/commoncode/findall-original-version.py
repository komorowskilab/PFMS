#
#  findall.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    pass

from commoncode import *
import sys, math

versionString = "%s: version 3.141592653589" % sys.argv[0]
print versionString

maxSpacing = 50
minHits = 4
minRatio = 4
minPlusRatio = 0.25
maxPlusRatio = 0.75
leftPlusRatio = 0.3
minPeak = 0.5
cachePages = -1
doPvalue = True
pValueType = 'self'
trimString = '10%'
shiftValue = 0
reportshift = False

if len(sys.argv) < 5:
    print "\nusage: python $ERANGEPATH/%s label samplerdsfile regionoutfile"  % sys.argv[0]
    print "[-control controlrdsfile] [-minimum minHits] [-ratio minRatio]"
    print "[-spacing maxSpacing] [-listPeak] [-shift #bp | learn] [-learnFold num]"  
    print "[-noshift] [-autoshift] [-reportshift] [-nomulti] [-minPlus fraction]"
    print "[-maxPlus fraction] [-leftPlus fraction] [-minPeak RPM] [-raw]"
    print "[-revbackground] [-pvalue self|back|none] [-nodirectionality]"
    print "[-strandfilter plus/minus] [-trimvalue percent] [-notrim]"
    print "[-cache pages] [-log altlogfile] [-flag aflag] [-append] [-RNA]\n"
    print "\twhere values in brackets are optional and label is an arbitrary string.\n"
    print "\t\tUse -ratio (default %d fold) to set the minimum fold enrichment" % minRatio
    print "\tover the control, -minimum (default %d) is the minimum number of reads" % minHits
    print "\t(RPM) within the region, and -spacing (default readlen) to set the maximum" 
    print "\tdistance between reads in the region. -listPeak lists the peak of the "
    print "\tregion. Peaks mut be higher than -minPeak (default %.1f RPM).\n" % minPeak
    print "\t\tPvalues are calculated from the sample (change with -pvalue),"
    print "\tunless the -revbackground flag and a control RDS file are provided.\n"
    print "\t\tBy default, all numbers and parameters are on a reads per "
    print "\tmillion (RPM) basis. -raw will treat all settings, ratios and reported"
    print "\t numbers as raw counts rather than RPM. Use -notrim to turn off region" 
    print "\ttrimming and -trimvalue to control trimming (default %s of peak signal)\n" % trimString
    print "\t\tThe peak finder uses minimal directionality information that can"
    print "\tbe turned off with -nodirectionality ; the fraction of + strand reads "
    print "\trequired to be to the left of the peak (default %.2f) can be set with " % leftPlusRatio
    print "\t-leftPlus ; -minPlus and -maxPlus change the minimum/maximum fraction"
    print "\tof plus reads in a region, which (defaults %.2f and %.2f, respectively).\n" % (minPlusRatio, maxPlusRatio)
    print "\t\tUse -shift to shift reads either by half the expected "
    print "\tfragment length (default 0 bp) or '-shift learn ' to learn the shift "
    print "\tbased on the first chromosome. If you prefer to learn the shift " 
    print "\tmanually, use -autoshift to calculate a per-region shift value, which "
    print "\tscan be reported using -reportshift. -strandfilter should only be used"
    print "\twhen explicitely calling unshifted stranded peaks from non-ChIP-seq " 
    print "\tdata such as directional RNA-seq. regionoutfile is written over by " 
    print "\tdefault unless given the -append flag.\n" 
    sys.exit(1)

factor = sys.argv[1]
hitfile = sys.argv[2]
outfilename = sys.argv[3]
logfilename = 'findall.log'
doAppend = False
stranded = 'both'

if '-append' in sys.argv:
    doAppend = True

if '-minimum' in sys.argv:
    minHits = float(sys.argv[sys.argv.index('-minimum') + 1])

if '-ratio' in sys.argv:
    minRatio = float(sys.argv[sys.argv.index('-ratio') + 1])

if '-spacing' in sys.argv:
    maxSpacing = int(sys.argv[sys.argv.index('-spacing') + 1])

if '-autoshift' in sys.argv:
    shiftValue = 'auto'
    shiftDict = {}

if '-shift' in sys.argv:
    try:
        shiftValue = int(sys.argv[sys.argv.index('-shift') + 1])
    except:
        try:
            if sys.argv[sys.argv.index('-shift') + 1] == 'learn':
                shiftValue = 'learn'
                print "Will try to learn shift"
        except:
            pass

if '-noshift' in sys.argv:
    shiftValue = 0
    
if '-reportshift' in sys.argv:
    reportshift = True

listPeak = False
if '-listPeak' in sys.argv:
    listPeak = True

normalize = True
if '-raw' in sys.argv:
    normalize = False
else:
    print "Normalizing to RPM"

if '-minPeak' in sys.argv:
    minPeak = float(sys.argv[sys.argv.index('-minPeak') + 1])

if '-minPlus' in sys.argv:
    minPlusRatio = float(sys.argv[sys.argv.index('-minPlus') + 1])

if '-maxPlus' in sys.argv:
    maxPlusRatio = float(sys.argv[sys.argv.index('-maxPlus') + 1])

doDirectionality = True
if '-nodirectionality' in sys.argv:
    doDirectionality = False

trimValue = 0.1
if '-trimvalue' in sys.argv:
    trimValue = float(sys.argv[sys.argv.index('-trimvalue') + 1]) / 100.
    trimString = '%2.1f' % (100. * trimValue)
    trimString += '\%'
    
doTrim = True
if '-notrim' in sys.argv:
    trimString = 'none'
    doTrim = False
    
if '-leftPlus' in sys.argv:
    leftPlusRatio = float(sys.argv[sys.argv.index('-leftPlus') + 1])

doRevBackground = False
if '-revbackground' in sys.argv:
    print "Swapping IP and background to calculate FDR"
    doRevBackground = True
    pValueType = 'back'

doControl = False
if '-control' in sys.argv:
    mockfile = sys.argv[sys.argv.index('-control') + 1]
    doControl = True

if '-pvalue' in sys.argv:
    ptype = sys.argv[sys.argv.index('-pvalue') + 1].upper()
    if ptype == 'NONE':
        doPvalue = False
        pValueType = 'none'
    elif ptype == 'SELF':
        pValueType = 'self'
    elif ptype == 'BACK':
        if doControl and doRevBackground:
            pValueType = 'back'
        else:
            print "must have a control dataset and -revbackground for pValue type 'back'"
    else:
        print "could not use pValue type : %s" % ptype

doCache = False
if '-cache' in sys.argv:
    doCache = True
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])

withFlag = ''
if '-flag' in sys.argv:
    withFlag = sys.argv[sys.argv.index('-flag') + 1]
    print "restrict to flag = %s" % withFlag

noMulti = False
useMulti = True
if '-nomulti' in sys.argv:
    print "using unique reads only"
    noMulti = True
    useMulti = False

if '-RNA' in sys.argv:
    print "using settings appropriate for RNA: -nodirectionality -notrim -noshift"
    shiftValue = 0
    doTrim = False
    doDirectionality = False
    
if '-strandfilter' in sys.argv:
    nextArg = sys.argv[sys.argv.index('-strandfilter') + 1]
    if nextArg == 'plus':
        stranded = '+'
        minPlusRatio = 0.9
        maxPlusRatio = 1.0
        print "only analyzing reads on the plus strand"
    elif nextArg == 'minus':
        stranded = '-'
        minPlusRatio = 0.0
        maxPlusRatio = 0.1
        print "only analyzing reads on the minus strand"

stringency = 4
if '-learnFold' in sys.argv:
    nextArg = sys.argv[sys.argv.index('-learnFold') + 1]
    if float(nextArg) > 1.0:
        stringency = float(nextArg)
        
if '-log' in sys.argv:
    logfilename = sys.argv[sys.argv.index('-log') + 1]

hitRDS = readDataset(hitfile, verbose = False)
readlen = hitRDS.getReadSize()
if '-RNA' in sys.argv:
    maxSpacing = readlen

writeLog(logfilename, versionString, string.join(sys.argv[1:]))

print "\nenforceDirectionality=%s listPeak=%s nomulti=%s cache=%s " % (doDirectionality, listPeak, noMulti, doCache)
print "spacing<%d minimum>%.1f ratio>%.1f minPeak=%.1f\ttrimmed=%s\tstrand=%s" % (maxSpacing, minHits, minRatio, minPeak, trimString, stranded)
try:
    print "minPlus=%.2f maxPlus=%.2f leftPlus=%.2f shift=%d pvalue=%s" % (minPlusRatio, maxPlusRatio, leftPlusRatio, shiftValue, pValueType)
except:
    print "minPlus=%.2f maxPlus=%.2f leftPlus=%.2f shift=%s pvalue=%s" % (minPlusRatio, maxPlusRatio, leftPlusRatio, shiftValue, pValueType)

if doControl:
    print "\ncontrol:" 
    mockRDS = readDataset(mockfile, verbose = True, cache=doCache)
    #sqlite default_cache_size is 2000 pages
    if cachePages > mockRDS.getDefaultCacheSize():
        mockRDS.setDBcache(cachePages)

print "\nsample:" 
hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
#sqlite default_cache_size is 2000 pages
if cachePages > hitRDS.getDefaultCacheSize():
    hitRDS.setDBcache(cachePages)
print 

hitDictSize = 1
hitRDSsize = len(hitRDS) / 1000000.
if doControl:
    mockDictSize = 1
    mockRDSsize = len(mockRDS) / 1000000.

if normalize:
    if doControl:
        mockSampleSize = mockRDSsize
    hitSampleSize = hitRDSsize

if doAppend:
    outfile = open(outfilename, 'a')
else:
    outfile = open(outfilename, 'w')

outfile.write('#ERANGE %s\n' % versionString)
if doControl:
    outfile.write('#enriched sample:\t%s (%.1f M reads)\n#control sample:\t%s (%.1f M reads)\n' % (hitfile, hitRDSsize, mockfile, mockRDSsize))
else:
    outfile.write('#enriched sample:\t%s (%.1f M reads)\n#control sample: none\n' % (hitfile, hitRDSsize))
if withFlag != '':
    outfile.write('#restrict to Flag = %s\n' % withFlag)
outfile.write('#enforceDirectionality=%s listPeak=%s nomulti=%s cache=%s\n' % (doDirectionality, listPeak, noMulti, doCache))
outfile.write('#spacing<%d minimum>%.1f ratio>%.1f minPeak=%.1f trimmed=%s strand=%s\n' % (maxSpacing, minHits, minRatio, minPeak, trimString, stranded))
try:
    outfile.write('#minPlus=%.2f maxPlus=%.2f leftPlus=%.2f shift=%d pvalue=%s\n' % (minPlusRatio, maxPlusRatio, leftPlusRatio, shiftValue, pValueType))
except:
    outfile.write('#minPlus=%.2f maxPlus=%.2f leftPlus=%.2f shift=%s pvalue=%s\n' % (minPlusRatio, maxPlusRatio, leftPlusRatio, shiftValue, pValueType))
if normalize:
    isNormalized = 'RPM'
else:
    isNormalized = 'COUNT'
headline = '#regionID\tchrom\tstart\tstop\t%s\tfold\tmulti' % isNormalized
headline += '%'
if doDirectionality:
    headline += '\tplus%\tleftPlus%'
if listPeak:
    headline += '\tpeakPos\tpeakHeight'
if reportshift:
    headline += '\treadShift'
if doPvalue:
    headline += '\tpValue'
#print headline
outfile.write('%s\n' %  headline)

index = 0
total = 0

mIndex = 0
mTotal = 0

if minRatio < minPeak:
    minPeak = minRatio

hitChromList = hitRDS.getChromosomes()
if doControl:
    mockChromList = mockRDS.getChromosomes()
hitChromList.sort()

failedCounter = 0
for achrom in hitChromList:
    outregions = []
    allregions = []
    if achrom == 'chrM':
        continue
    if doControl and (achrom not in mockChromList):
        continue
    print "chromosome %s" % (achrom)
    previousHit = - 1 * maxSpacing
    currentHitList = [-1]
    currentWeightList = [0]
    currentReadList = []
    mockStart = 0
    numStarts = 0
    if withFlag == '':
        hitDict = hitRDS.getReadsDict(fullChrom=True, chrom=achrom, withWeight=True, doMulti=useMulti, findallOptimize=True, strand=stranded)
    else:
        hitDict = hitRDS.getReadsDict(fullChrom=True, chrom=achrom, flag=withFlag, withWeight=True, doMulti=useMulti, findallOptimize=True, strand=stranded)
    #mockRDS.memSync(chrom=achrom, index=True)
    maxCoord = hitRDS.getMaxCoordinate(achrom, doMulti=useMulti)
    chromIndex = 1
    if shiftValue == 'learn':
        print "learning shift.... will need at least 30 training sites"
        learnPreviousHit = -1 * maxSpacing
        learnHitList = [-1]
        learnWeightList = [0]
        learnReadList = []
        learnDict = {}
        learnCount = 0
        for (pos, sense, weight) in hitDict[achrom]:
            if abs(pos - learnPreviousHit) > maxSpacing or pos == maxCoord:
                sumAll = sum(learnWeightList)
                if normalize:
                    sumAll /= hitSampleSize
                regionStart = learnHitList[0]
                regionStop = learnHitList[-1]
                # we're going to require stringent settings
                if sumAll >= stringency * minHits and numStarts > stringency * minRatio and (regionStop - regionStart) > stringency * readlen:
                    if doControl:
                        numMock = 1. + mockRDS.getCounts(achrom, regionStart, regionStop, uniqs=True, multi=useMulti, splices=False, reportCombined=True)
                        if normalize:
                            numMock /= mockSampleSize
                        foldRatio = sumAll / numMock
                    else:
                        foldRatio = minRatio
                    if foldRatio >= minRatio:
                        (topPos, numHits, smoothArray, numPlus, localshift) = findPeak(learnReadList, regionStart, regionStop - regionStart, readlen, doWeight=True, shift='auto', returnShift=True)
                        try:
                            learnDict[localshift] += 1
                        except:
                            learnDict[localshift] = 1
                        learnCount += 1
                learnHitList = []
                learnWeightList = []
                learnReadList = []
                learnStarts = 0
            if pos not in currentHitList:
                numStarts += 1
            learnHitList.append(pos)
            learnWeightList.append(weight)
            learnReadList.append((pos, sense, weight))
            learnPreviousHit = pos
        bestShift = 0
        bestCount = 0
        outline = "#learn: stringency=%.2f min_signal=%2.f min_ratio=%.2f min_region_size=%d\n#number of training examples: %d" % (stringency,stringency * minHits,stringency * minRatio,stringency * readlen,learnCount)
        print outline
        writeLog(logfilename, versionString, outfilename + outline)
        if learnCount < 30:
            outline = "#too few training examples to pick a shiftValue - defaulting to 0\n"
            outline = "#consider picking a lower minimum or threshold"
            print outline
            writeLog(logfilename, versionString, outfilename + outline)
            shiftValue = 0
        else:
            for shift in sorted(learnDict):
                if learnDict[shift] > bestCount:
                    bestShift = shift
                    bestCount = learnDict[shift]
            shiftValue = bestShift
            print learnDict
        outline = "#picked shiftValue to be %d" % shiftValue
        print outline
        outfile.write(outline + '\n')
        writeLog(logfilename, versionString, outfilename + outline)
    
    for (pos, sense, weight) in hitDict[achrom]:
        chromIndex += 1
        if abs(pos - previousHit) > maxSpacing or pos == maxCoord:
            sumAll = sum(currentWeightList)
            if normalize:
                sumAll /= hitSampleSize
            regionStart = currentHitList[0]
            regionStop = currentHitList[-1]
            allregions.append(int(sumAll))
            if sumAll >= minHits and numStarts > minRatio and (regionStop - regionStart) > readlen:
                sumMulti = 0.
                if doControl:
                    numMock = 1. + mockRDS.getCounts(achrom, regionStart, regionStop, uniqs=True, multi=useMulti, splices=False, reportCombined=True)
                    if normalize:
                        numMock /= mockSampleSize
                    foldRatio = sumAll / numMock
                else:
                    foldRatio = minRatio
                if foldRatio >= minRatio:
                    # first pass, with absolute numbers
                    if doDirectionality:
                        (topPos, numHits, smoothArray, numPlus, numLeft, localshift) = findPeak(currentReadList, regionStart, regionStop - regionStart, readlen, doWeight=True, leftPlus=doDirectionality, shift=shiftValue, returnShift=True)
                    else:
                        (topPos, numHits, smoothArray, numPlus, localshift) = findPeak(currentReadList, regionStart, regionStop - regionStart, readlen, doWeight=True, shift=shiftValue, returnShift=True)
                    bestPos = topPos[0]
                    peakScore = smoothArray[bestPos]
                    if normalize:
                        peakScore /= hitSampleSize
                    if doTrim:
                        minSignalThresh = trimValue * peakScore 
                        start = 0
                        stop = regionStop - regionStart - 1
                        startFound = False
                        while not startFound:
                            if smoothArray[start] >= minSignalThresh or start == bestPos:
                                startFound = True
                            else:
                                start += 1
                        stopFound = False
                        while not stopFound:
                            if smoothArray[stop] >= minSignalThresh or stop == bestPos:
                                stopFound = True
                            else:
                                stop -= 1
                        regionStop = regionStart + stop
                        regionStart += start
                        try:
                            if doDirectionality:
                                (topPos, sumAll, smoothArray, numPlus, numLeft) = findPeak(currentReadList, regionStart, regionStop - regionStart, readlen, doWeight=True, leftPlus=doDirectionality, shift=localshift)
                            else:
                                (topPos, sumAll, smoothArray, numPlus) = findPeak(currentReadList, regionStart, regionStop - regionStart, readlen, doWeight=True, shift=localshift)
                        except:
                            continue
                        if normalize:
                            sumAll /= hitSampleSize
                        if doControl:
                            numMock = 1. + mockRDS.getCounts(achrom, regionStart, regionStop, uniqs=True, multi=useMulti, splices=False, reportCombined=True)
                            if normalize:
                                numMock /= mockSampleSize
                            foldRatio = sumAll / numMock
                        else:
                            foldRation = minRatio
                        sumMulti = hitRDS.getCounts(achrom, regionStart, regionStop, uniqs=False, multi=useMulti, splices=False, reportCombined=True)
                        # just in case it changed, use latest data
                        try:
                            bestPos = topPos[0]
                            peakScore = smoothArray[bestPos]
                            if normalize:
                                peakScore /= hitSampleSize
                        except:
                            continue
                    else:
                        sumMulti = sum(currentWeightList) - currentWeightList.count(1.0)
                    # normalize to RPM
                    if normalize:
                        sumMulti /= hitSampleSize
                    try:
                        multiP = 100. * (sumMulti / sumAll)
                    except:
                        break
                    if noMulti:
                        multiP = 0.
                    # check that we still pass threshold
                    if sumAll >= minHits and  foldRatio >= minRatio and (regionStop - regionStart) > readlen:
                        plusRatio = float(numPlus)/numHits
                        #print ('considering', achrom, regionStart, regionStop, numHits, numPlus, numLeft, plusRatio)
                        if peakScore >= minPeak and minPlusRatio <= plusRatio <= maxPlusRatio:
                            peak = ''
                            if listPeak:
                                peak= '\t%d\t%.1f' % (regionStart + bestPos, peakScore)
                            if doDirectionality:
                                if leftPlusRatio < numLeft / numPlus:
                                    index += 1
                                    plusP = plusRatio * 100.
                                    leftP = 100. * numLeft / numPlus
                                    # we have a region that passes all criteria
                                    outregions.append((factor, index, achrom, regionStart, regionStop + readlen - 1, sumAll, foldRatio, multiP, plusP, leftP, peak, localshift))
                                    #outline = '%s%d\t%s\t%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f%s' % (factor, index, achrom, regionStart, regionStop + readlen - 1, sumAll, foldRatio, multiP, plusP, leftP, peak)
                                    #print outline
                                    #outfile.write('%s\n' %  outline)
                                    total += sumAll
                                else:
                                    #print ('failed directionality', achrom, regionStart, regionStop, numHits, numPlus, numLeft)
                                    failedCounter += 1
                            else:
                                index += 1
                                # we have a region, but didn't check for directionality
                                outregions.append((factor, index, achrom, regionStart, regionStop + readlen - 1, sumAll, foldRatio, multiP, peak, localshift))
                                #outline = '%s%d\t%s\t%d\t%d\t%.1f\t%.1f\t%.1f%s' % (factor, index, achrom, regionStart, regionStop + readlen - 1, sumAll, foldRatio, multiP, peak)
                                #print outline
                                #outfile.write('%s\n' %  outline)
                                total += sumAll
            currentHitList = []
            currentWeightList = []
            currentReadList = []
            numStarts = 0
        if pos not in currentHitList:
            numStarts += 1
        currentHitList.append(pos)
        currentWeightList.append(weight)
        currentReadList.append((pos, sense, weight))
        previousHit = pos
    
    if not doRevBackground:
        if doPvalue:
            allregions.sort()
            regionallSize = float(len(allregions))
            try:
                poissonmean = sum(allregions) / regionallSize
            except:
                poissonmean = 0
            print "Poisson n=%d, p=%f" % (regionallSize, poissonmean)
            p = math.exp(-poissonmean)
        print headline
        for region in outregions:
            # iterative poisson from http://stackoverflow.com/questions/280797?sort=newest
            if doPvalue:
                pValue = p
                sumAll = int(region[5])
                for i in xrange(sumAll):
                    pValue *= poissonmean
                    pValue /= i+1
            try:
                if reportshift:
                    outline = '%s%d\t%s\t%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f%s\t%d' % region
                else:
                    outline = '%s%d\t%s\t%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f%s' % region[:-1]
            except:
                if reportshift:
                    outline = '%s%d\t%s\t%d\t%d\t%.1f\t%.1f\t%.1f%s\t%d' % region
                else:
                    outline = '%s%d\t%s\t%d\t%d\t%.1f\t%.1f\t%.1f%s' % region[:-1]
            if doPvalue:
                outline += '\t%1.2g' % pValue
            print outline
            outfile.write(outline + '\n')
        continue
    #now do background swapping the two samples around
    print "calculating background..."
    previousHit = - 1 * maxSpacing
    currentHitList = [-1]
    currentWeightList = [0]
    currentReadList = []
    backregions = []
    mockStart = 0
    numStarts = 0
    hitDict = mockRDS.getReadsDict(fullChrom=True, chrom=achrom, withWeight=True, doMulti=useMulti, findallOptimize=True)
    #mockRDS.memSync(chrom=achrom, index=True)
    maxCoord = mockRDS.getMaxCoordinate(achrom, doMulti=useMulti)
    for (pos, sense, weight) in hitDict[achrom]:
        if abs(pos - previousHit) > maxSpacing or pos == maxCoord:
            sumAll = sum(currentWeightList)
            if normalize:
                sumAll /= mockSampleSize
            regionStart = currentHitList[0]
            regionStop = currentHitList[-1]
            backregions.append(int(sumAll))
            if sumAll >= minHits and numStarts > minRatio and (regionStop - regionStart) > readlen:
                numMock = 1. + hitRDS.getCounts(achrom, regionStart, regionStop, uniqs=True, multi=useMulti, splices=False, reportCombined=True)
                if normalize:
                    numMock /= hitSampleSize
                foldRatio = sumAll / numMock
                if foldRatio >= minRatio:
                    # first pass, with absolute numbers
                    if doDirectionality:
                        (topPos, numHits, smoothArray, numPlus, numLeft, backshift) = findPeak(currentReadList, regionStart, regionStop - regionStart, readlen, doWeight=True, leftPlus=doDirectionality, shift=shiftValue, returnShift=True)
                    else:
                        (topPos, numHits, smoothArray, numPlus, backshift) = findPeak(currentReadList, regionStart, regionStop - regionStart, readlen, doWeight=True, shift=shiftValue, returnShift=True)
                    bestPos = topPos[0]
                    peakScore = smoothArray[bestPos]
                    if normalize:
                        peakScore /= mockSampleSize
                    if doTrim:
                        minSignalThresh = peakScore / 20.
                        start = 0
                        stop = regionStop - regionStart - 1
                        startFound = False
                        while not startFound:
                            if smoothArray[start] >= minSignalThresh or start == bestPos:
                                startFound = True
                            else:
                                start += 1
                        stopFound = False
                        while not stopFound:
                            if smoothArray[stop] >= minSignalThresh or stop == bestPos:
                                stopFound = True
                            else:
                                stop -= 1
                        regionStop = regionStart + stop
                        regionStart += start
                        try:
                            if doDirectionality:
                                (topPos, sumAll, smoothArray, numPlus, numLeft) = findPeak(currentReadList, regionStart, regionStop - regionStart, readlen, doWeight=True, leftPlus=doDirectionality, shift=backshift)
                            else:
                                (topPos, sumAll, smoothArray, numPlus) = findPeak(currentReadList, regionStart, regionStop - regionStart, readlen, doWeight=True, shift=backshift)
                        except:
                            continue
                        numMock = 1. + hitRDS.getCounts(achrom, regionStart, regionStop, uniqs=True, multi=useMulti, splices=False, reportCombined=True)
                        if normalize:
                            sumAll /= mockSampleSize
                            numMock /= hitSampleSize
                        foldRatio = sumAll / numMock
                        # just in case it changed, use latest data
                        try:
                            bestPos = topPos[0]
                            peakScore = smoothArray[bestPos]
                        except:
                            continue
                        # normalize to RPM
                        if normalize:
                            peakScore /= mockSampleSize
                    # check that we still pass threshold
                    if sumAll >= minHits and  foldRatio >= minRatio and (regionStop - regionStart) > readlen:
                        plusRatio = float(numPlus)/numHits
                        #print ('considering', achrom, regionStart, regionStop, numHits, numPlus, numLeft, plusRatio)
                        if peakScore >= minPeak and minPlusRatio <= plusRatio <= maxPlusRatio:
                            if doDirectionality:
                                if leftPlusRatio < numLeft / numPlus:
                                    mIndex += 1
                                    mTotal += sumAll
                                else:
                                    #print ('failed directionality', achrom, regionStart, regionStop, numHits, numPlus, numLeft)
                                    failedCounter += 1
                            else:
                                # we have a region, but didn't check for directionality
                                mIndex += 1
                                mTotal += sumAll
            currentHitList = []
            currentWeightList = []
            currentReadList = []
            numStarts = 0
        if pos not in currentHitList:
            numStarts += 1
        currentHitList.append(pos)
        currentWeightList.append(weight)
        currentReadList.append((pos, sense, weight))
        previousHit = pos
    print mIndex, mTotal
    if doPvalue:
        if pValueType == 'self':
            allregions.sort()
            regionallSize = float(len(allregions))
            try:
                poissonmean = sum(allregions) / regionallSize
            except:
                poissonmean = 0
            print "Poisson n=%d, p=%f" % (regionallSize, poissonmean)
        else:
            backregions.sort()
            backregionSize = float(len(backregions))
            try:
                poissonmean = sum(backregions) / backregionSize
            except:
                poissonmean = 0
            print "Poisson n=%d, p=%f" % (backregionSize, poissonmean)        
        p = math.exp(-poissonmean)
    print headline
    for region in outregions:
        # iterative poisson from http://stackoverflow.com/questions/280797?sort=newest
        if doPvalue:
            pValue = p
            sumAll = int(region[5])
            for i in xrange(sumAll):
                pValue *= poissonmean
                pValue /= i+1
        if shiftValue == 'auto' and reportshift:
            try:
                shiftDict[region[-1]] += 1
            except:
                shiftDict[region[-1]] = 1
        try:
            if reportshift:
                outline = '%s%d\t%s\t%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f%s\t%d' % region
            else:
                outline = '%s%d\t%s\t%d\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f%s' % region[:-1]
        except:
            if reportshift:
                outline = '%s%d\t%s\t%d\t%d\t%.1f\t%.1f\t%.1f%s\t%d' % region
            else:
                outline = '%s%d\t%s\t%d\t%d\t%.1f\t%.1f\t%.1f%s' % region[:-1]
        if doPvalue:
            outline += '\t%1.2g' % pValue
        print outline
        outfile.write(outline + '\n')
    
footer = '#stats:\t%.1f RPM in %d regions\n' % (total, index)
if doDirectionality:
    footer += '#\t\t%d additional regions failed directionality filter\n' % failedCounter

if doRevBackground:
    try:
        percent = 100. * (float(mIndex)/index)
    except:
        percent = 0.
    if percent > 100.:
        percent = 100.
    footer += '#%d regions (%.1f RPM) found in background (FDR = %.2f percent)\n' % (mIndex, mTotal, percent)
if shiftValue == 'auto' and reportshift:
    bestShift = 0
    bestCount = 0
    for shift in sorted(shiftDict):
        if shiftDict[shift] > bestCount:
            bestShift = shift
            bestCount = shiftDict[shift]
    footer += '#mode of shift values: %d\n' % bestShift

print footer
outfile.write(footer)
outfile.close()

writeLog(logfilename, versionString, outfilename + footer.replace('\n#',' | ')[:-1])
