#
#  getSNPs.py
#  ENRAGE
#
# Get the matches and mismatches from the RDS file, and calulate the SNP thresholds uniqStartMin (Sl * readlength) and and totalRatio (Cl). 
# For each mismatch, choose the base change that occur most frequently (ie: has the highest number
# of indpendent reads)
# Threshold of Sl and Cl are from user input
# Sl = # of independent reads supporting a base change at position S 
# Cl = total # of all reads suppoting a base change at position S / # of all # reads that pass through position S
# Originally written by: Wendy Lee
# Last modified: May 11th, 2009 by Ali Mortazavi

import sys
from commoncode import *

print '%s: version 3.5' % sys.argv[0]

try:
    import psyco
    psyco.full()
except:
    print "psyco is not running"
    pass

def getMatchDict(rds, achrom, withSplices = True):
    spliceDict = {}
    newDict = {}
    finalDict = {}
    try:
        newDict = rds.getReadsDict(fullChrom=True, bothEnds=True, noSense=True, chrom=achrom)
    except:
        newDict[achrom] = []
    for (start, stop) in newDict[achrom]:
        try:
            finalDict[start].append(stop)
        except:
            finalDict[start] = [stop]
            
    if withSplices:
        try:
            spliceDict = rds.getSplicesDict(noSense=True, fullChrom=True, chrom=achrom, splitRead=True)
        except:
            spliceDict[achrom] = []
        for (start, stop) in spliceDict[achrom]:
            try:
                finalDict[start].append(stop)
            except:
                finalDict[start] = [stop]
    
    return finalDict

def getMismatchDict(theRDS, chrom_in, withSplices = True):
    mainDict = {} 
    spDict = theRDS.getMismatches(mischrom=chrom_in, useSplices=withSplices)
    for (start, change_at, change_base, change_from) in spDict[chrom_in]:
        change = change_base + "-" + change_from
        try:
            (ucount, totalcount, back, uniqBaseDict, totalBaseDict) = mainDict[change_at]
            pos = str(start) + ':' + change
            try: 
                totalBaseDict[change] += 1
            except: 
                totalBaseDict[change] = 1
            if pos not in back:
                ucount += 1 # unique read count
                try:
                    uniqBaseDict[change] += 1 # dict contains total unique read counts
                except:
                    uniqBaseDict[change] = 1
                mainDict[change_at] = [ucount, totalcount + 1, back + "," + pos, uniqBaseDict, totalBaseDict]
            else: # if this is not a unique read
                mainDict[change_at] = [ucount, totalcount + 1, back, uniqBaseDict, totalBaseDict]
        except: 
            mainDict[change_at] = [1,1, str(start) + ":" + change, {change:1}, {change:1}] 

    return mainDict

#Main program 
if len(sys.argv) < 5:
    print 'usage: python %s samplerdsfile uniqStartMin totalRatioMin outfile [-nosplices] [-enforceChr] [-cache pages]' % sys.argv[0]
    print 'where\n\tuniqStartMin = # of independent reads supporting a base change at position S'
    print '\ttotalRatioMin = total # of reads suppoting a base change at position S / total # reads that pass through position S'
    sys.exit(1)

outStr = "" 
hitfile = sys.argv[1]
uniqStartMin = float(sys.argv[2])
totalRatioMin = float(sys.argv[3])
writeLog('snp.log', sys.argv[0], "rdsfile: %s uniqStartMin: %1.2f totalRatioMin: %1.2f" % (hitfile, uniqStartMin, totalRatioMin))

outfilename = sys.argv[4]
cachePages = 0
doCache = False
if '-cache' in sys.argv:
    cachePages = int(sys.argv[sys.argv.index('-cache') + 1])
if cachePages > 0:
    doCache = True

forceChr = False
if '-enforceChr' in sys.argv:
    forceChr = True

doSplices = True
if '-nosplices' in sys.argv:
    doSplices = False

hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
if cachePages > 20000:
    hitRDS.setDBcache(cachePages)
readLength = hitRDS.getReadSize() 
outfile  = open(outfilename, 'w')
chromList = hitRDS.getChromosomes()

snpCount = 0
header = "#Sl\tCl\tchrom\tpos\tmatch\tuniqMis\t\ttotalMis\tchange" 
outfile.write(header + '\n')
for chrom_in in chromList:
    if forceChr:
        if chrom_in[:3] != 'chr':
            continue
    matchDict = getMatchDict(hitRDS, chrom_in, doSplices)
    print 'got match dict for %s ' % chrom_in
    mismatchDict = getMismatchDict(hitRDS, chrom_in, doSplices)
    print 'got mismatch dict for %s ' % chrom_in
    start = 0
    misPositions = mismatchDict.keys()
    misPositions.sort()
    print header
    for position in misPositions:
        (uniqCnt, totalCount, back, uniqBaseDict, totalBaseDict) = mismatchDict[position]
        highestCount = 0
        highestBaseChange = 'N-N'
        highestTotalCount = 0
        for baseChange in uniqBaseDict:
            if totalBaseDict[baseChange] > highestTotalCount:
                highestBaseChange = baseChange
                highestCount = uniqBaseDict[baseChange]
                highestTotalCount = totalBaseDict[baseChange]
        Cl = 0.
        matchCount = 0
        if highestCount >= uniqStartMin:
            #(skip, matchCount) = checkBetween(matchList[start:], position) 
            for matchpos in xrange(position - readLength + 1, position + 1):
                try:
                    matchCount += len([mstop for mstop in matchDict[matchpos] if position <= mstop])
                except:
                    pass
            matchCount -= totalCount
            if matchCount < 0:
                matchCount = 0
            Sl = highestCount/float(readLength)
            Cl = highestTotalCount/float(highestTotalCount + matchCount)
            if Cl >= totalRatioMin:
                outline = "%1.2f\t%1.2f\t%s\t%d\t%d\t%d\t\t%d\t%s" % (Sl, Cl, chrom_in, position, matchCount, highestCount, highestTotalCount, highestBaseChange)
                print outline
                outfile.write(outline + "\n")
                outfile.flush() 
                snpCount += 1
outfile.close()

writeLog('snp.log', sys.argv[0], "%d candidate SNPs\n" % snpCount)

