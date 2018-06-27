#
#  weightMultireads.py
#  ENRAGE
#

#  Created by Ali Mortazavi on 10/02/08.
#

try:
    import psyco
    psyco.full()
except:
    pass

from commoncode import *
import sys, time

print '%s: version 3.1' % sys.argv[0]

if len(sys.argv) < 2:
    print 'usage: python %s rdsfile [-radius bp] [-noradius] [-usePairs maxDist] [-verbose] [-cache pages]' % sys.argv[0]
    sys.exit(1)

rdsfile = sys.argv[1]

radius = 100
doRadius = True
if '-noradius' in sys.argv:
    doRadius = False

if '-radius' in sys.argv:
    radius = int(sys.argv[sys.argv.index('-radius')+1])
    doRadius = True

usePairs = False
pairDist = 200
if '-usePairs' in sys.argv:
    usePairs = True
    try:
        pairDist = int(sys.argv[sys.argv.index('-usePairs')+1])
    except:
        pairDist = 200

tooFar = pairDist * 10

verbose = False
if '-verbose' in sys.argv:
    verbose = True
    
cachePages = 1
doCache = False
if '-cache' in sys.argv:
    doCache = True
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])

RDS = readDataset(rdsfile, verbose = True, cache=doCache)
readlen = RDS.getReadSize()
halfreadlen = readlen / 2

if cachePages > RDS.getDefaultCacheSize():
    RDS.setDBcache(cachePages)

if verbose:
    print time.ctime()
multiIDs = RDS.getReadIDs(uniqs=False,multi=True)
if verbose:
    print 'got multiIDs ', time.ctime()

fixedPair = 0
fixedReads = []
if usePairs:
    print 'doing pairs with pairDist = %d' % pairDist
    uidDict = {}
    midDict = {}
    jointList = []
    bothMultiList = []
    mainIDList = []
    guDict = {}
    muDict = {}
    
    if RDS.dataType == 'RNA':
        uniqIDs = RDS.getReadIDs(uniqs=True,multi=False,splices=True)
    else:
        uniqIDs = RDS.getReadIDs(uniqs=True,multi=False,splices=False)
    
    if verbose:
        print 'got uniqIDs ', time.ctime()
    
    for readID in uniqIDs:
        (mainID, pairID) = readID.split('/')
        try:
            uidDict[mainID].append(pairID)
        except:
            uidDict[mainID] = [pairID]
            mainIDList.append(mainID)
    if verbose:
        print 'uidDict all' , len(uidDict), time.ctime()
    
    for mainID in mainIDList:
        if len(uidDict[mainID]) == 2:
            del uidDict[mainID]
    if verbose:
        print 'uidDict first candidates ', len(uidDict), time.ctime()
    
    for readID in multiIDs:
        (frontID, multiplicity) = readID.split('::')
        (mainID, pairID) = frontID.split('/')
        try:
            if pairID not in midDict[mainID]:
                midDict[mainID].append(pairID)
        except:
            midDict[mainID] = [pairID]
    if verbose:
        print 'all multis ', len(midDict), time.ctime()
    
    mainIDList = uidDict.keys()
    for mainID in mainIDList:
        if mainID not in midDict:
            del uidDict[mainID]
    if verbose:
        print 'uidDict actual candidates ', len(uidDict), time.ctime()
    
    for readID in midDict:
        listLen = len(midDict[readID])
        if listLen == 1:
            if readID in uidDict:
                jointList.append(readID)
        elif listLen == 2:
            bothMultiList.append(readID)
    if verbose:
        print 'joint ', len(jointList), time.ctime()
        print 'bothMulti ', len(bothMultiList), time.ctime()
    del uidDict
    del midDict
    del mainIDList
    del uniqIDs
    
    uniqDict = RDS.getReadsDict(noSense=True, withChrom=True, withPairID=True, doUniqs=True, readIDDict=True)
    if verbose:
        print 'got uniq dict', len(uniqDict), time.ctime()
    if RDS.dataType == 'RNA':
        spliceDict = RDS.getSplicesDict(noSense=True, withChrom=True, withPairID=True, readIDDict=True)
        if verbose:
            print 'got splice dict', len(spliceDict), time.ctime()
    
    for readID in jointList:
        try:
            guDict[readID] = uniqDict[readID][0]
        except:
            if RDS.dataType == 'RNA':
                guDict[readID] = spliceDict[readID][0]
    del uniqDict
    del spliceDict
    if verbose:
        print 'guDict actual ', len(guDict), time.ctime()
    
    multiDict = RDS.getReadsDict(noSense=True, withChrom=True, withPairID=True, doUniqs=False, doMulti=True, readIDDict=True)
    if verbose:
        print 'got multi dict', len(multiDict), time.ctime()
    
    for readID in jointList:
        muDict[readID] = multiDict[readID]
    
    for readID in bothMultiList:
        muDict[readID] = multiDict[readID]
    del multiDict
    if verbose:
        print 'muDict actual ', len(muDict), time.ctime()
    
    RDS.setSynchronousPragma('OFF')
    for readID in jointList:
        try:
            (ustart, uchrom, upair) = guDict[readID]
            ustop = ustart + readlen
        except:
            (ustart, lstop, rstart, ustop, uchrom, upair) = guDict[readID]
        muList = muDict[readID]
        muLen = len(muList)
        bestMatch = [tooFar] * muLen
        found = False
        for index in range(muLen):
            (mstart, mchrom, mpair) = muList[index]
            if uchrom != mchrom:
                continue
            if abs(mstart - ustart) < pairDist:
                bestMatch[index] = abs(mstart - ustart)
                found = True
            elif abs(mstart - ustop) < pairDist:
                bestMatch[index] = abs(mstart - ustop)
                found = True
        if found:
            theMatch = -1
            theDist = tooFar
            reweighList = []
            for index in range(muLen):
                if theDist > bestMatch[index]:
                    theMatch = index
                    theDist = bestMatch[index]
            theID = readID + '/' +  mpair
            for index in range(muLen):
                if index == theMatch:
                    score = 1 - (muLen - 1) / (100. * (muLen))
                else:
                    score = 1 / (100. * muLen)
                start = muList[index][0]
                chrom = 'chr' + muList[index][1]
                reweighList.append((round(score,3), chrom, start, theID))
            if theMatch > 0:
                RDS.reweighMultireads(reweighList)
                fixedPair += 1
                if verbose and fixedPair % 10000 == 1:
                    print 'fixed %d' % fixedPair
                    print guDict[readID]
                    print muDict[readID]
                    print reweighList
                fixedReads.append(theID)
    RDS.setSynchronousPragma('ON')
    
    del guDict
    del muDict
    print 'fixed %d pairs' % fixedPair
    print time.ctime()

skippedReads = 0
if doRadius:
    print 'doing uniq read radius with radius = %d' % radius
    multiDict = RDS.getReadsDict(noSense=True, withWeight=True, withChrom=True, withID=True, doUniqs=False, doMulti=True, readIDDict=True)
    print 'got multiDict'
    RDS.setSynchronousPragma('OFF')
    rindex = 0
    for readID in multiIDs:
        theID = readID
        if theID in fixedReads:
            skippedReads += 1
            continue
        if '::' in readID:
            (readID, multiplicity) = readID.split('::')
        #hitDict = RDS.getReadsDict(withWeight=True, withChrom=True, withID=True, doUniqs=False, doMulti=True, readIDDict=True, readLike=readID)
        scores = []
        coords = []
        for read in multiDict[readID]:
            (start, weight, rID, chrom) = read
            achrom = 'chr' + chrom
            regionStart = start + halfreadlen - radius
            regionStop = start + halfreadlen + radius 
            uniqs = RDS.getCounts(achrom, regionStart, regionStop, uniqs=True, multi=False, splices=False, reportCombined=True)
            scores.append(uniqs + 1)
            coords.append((achrom, start, theID))
        total = float(sum(scores))
        reweighList = []
        for index in range(len(scores)):
            reweighList.append((round(scores[index]/total,2), coords[index][0], coords[index][1], coords[index][2]))
        RDS.reweighMultireads(reweighList)
        rindex += 1
        if rindex % 10000 == 0:
            print rindex
    RDS.setSynchronousPragma('ON')
    if verbose:
        print 'skipped ', skippedReads
    print 'reweighted ', rindex

if doCache:
    RDS.saveCacheDB(rdsfile)

if verbose:
    print 'finished', time.ctime()
#for readID in hitDict:
#    reads = hitDict[readID]
    
