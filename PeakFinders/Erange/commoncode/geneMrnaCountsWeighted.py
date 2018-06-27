try:
    import psyco
    psyco.full()
except:
    print 'psyco not running'

import sys

print '%s: version 4.1' % sys.argv[0]

if len(sys.argv) < 5:
    print 'usage: python %s genome rdsfile uniqcountfile outfile [-stranded] [-uniq] [-multi] [-record] [-accept acceptfile] [-cache pages] [-verbose]' % sys.argv[0]
    sys.exit(1)

from commoncode import *
from cistematic.genomes import Genome
from cistematic.core.geneinfo import geneinfoDB
from cistematic.core import chooseDB, cacheGeneDB, uncacheGeneDB

genome = sys.argv[1]
hitfile =  sys.argv[2]
countfile = sys.argv[3]
outfilename = sys.argv[4]

withUniqs = False
if '-uniq' in sys.argv:
    withUniqs = True

withMulti = False
if '-multi' in sys.argv:
    withMulti = True

if (not withUniqs and not withMulti) or (withUniqs and withMulti):
    print "must have either one of -uniq or -multi set. Exiting"
    sys.exit(1)

recording = False
if '-record' in sys.argv:
    if withUniqs:
        print "ignoring -record with uniq reads"
    else:
        recording = True

ignoreSense = True
if '-stranded' in sys.argv:
    print "will track strandedness"
    ignoreSense = False
    
doVerbose = False
if '-verbose' in sys.argv:
    doVerbose = True
    
if '-accept' in sys.argv:
    acceptfile = sys.argv[sys.argv.index('-accept') + 1]
else:
    acceptfile = 'none'
    acceptDict = {}

if acceptfile != 'none':
    acceptDict = getMergedRegions(acceptfile, maxDist=0, keepLabel=True, verbose=True)

extendGenome = ''
replaceModels = False
if '-models' in sys.argv:
    extendGenome = sys.argv[sys.argv.index('-models') + 1]
    if '-replacemodels' in sys.argv:
        replaceModels = True
        print "will replace gene models with %s" % extendGenome
    else:
        print "will extend gene models with %s" % extendGenome

doCache = False
cachePages = 0
if '-cache' in sys.argv:
    cacheGeneDB(genome)
    hg = Genome(genome, dbFile=chooseDB(genome), inRAM=True)
    idb = geneinfoDB(cache=True)
    print '%s cached' % genome
    doCache = True
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])
else:
    hg = Genome(genome, inRAM=True)
    idb = geneinfoDB()

if extendGenome != '':
    hg.extendFeatures(extendGenome, replace = replaceModels)
    
hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
if cachePages > hitRDS.getDefaultCacheSize():
    hitRDS.setDBcache(cachePages)

readlen = hitRDS.getReadSize()

geneinfoDict = idb.getallGeneInfo(genome)
geneannotDict = hg.allAnnotInfo()
gidCount = {}
gidReadDict = {}

featuresByChromDict = getFeaturesByChromDict(hg, acceptDict)
gidList = hg.allGIDs()

gidList.sort()
for chrom in acceptDict:
    for (label, start, stop, length) in acceptDict[chrom]:
        if label not in gidList:
            gidList.append(label)

for gid in gidList:
    gidCount[gid] = 0
    gidReadDict[gid] = []

uniqueCountDict = {}
gidDict = {}
read2GidDict = {}

uniquecounts = open(countfile)
for line in uniquecounts:
    fields = line.strip().split()
    # add a pseudo-count here to ease calculations below
    uniqueCountDict[fields[0]] = float(fields[-1]) + 1
uniquecounts.close()

#infile = open(infilename)
outfile = open(outfilename,'w')

index = 0
if withMulti and not withUniqs:
    chromList = hitRDS.getChromosomes(table='multi', fullChrom=False)
else:
    chromList = hitRDS.getChromosomes(fullChrom=False)

for achrom in chromList:
    if achrom not in featuresByChromDict:
        continue
    print "\n" + achrom + " ",
    startFeature = 0
    fullchrom = 'chr' + achrom
    hitDict = hitRDS.getReadsDict(noSense=ignoreSense, fullChrom=True, chrom=fullchrom, withID=True, doUniqs=withUniqs, doMulti=withMulti)
    featList = featuresByChromDict[achrom]
    if ignoreSense:
        for (tagStart, tagReadID) in hitDict[fullchrom]:
            index += 1
            if index % 100000 == 0:
                print "read %d" % index,
            stopPoint = tagStart + readlen
            if startFeature < 0:
                startFeature = 0
            for (start, stop, gid, sense, ftype) in featList[startFeature:]:
                if tagStart > stop:
                    startFeature += 1
                    continue
                if start > stopPoint:
                    startFeature -= 100
                    break
                if start <= tagStart <= stop:
                    try:
                        gidReadDict[gid].append(tagReadID)
                        if tagReadID in read2GidDict:
                            if gid not in read2GidDict[tagReadID]:
                                read2GidDict[tagReadID].append(gid)
                        else:
                            read2GidDict[tagReadID] = [gid]
                        gidCount[gid] += 1
                    except:
                        print "gid %s not in gidReadDict" % gid
                    stopPoint = stop
    else:
        for (tagStart, tSense, tagReadID) in hitDict[fullchrom]:
            index += 1
            if index % 100000 == 0:
                print "read %d" % index,
            stopPoint = tagStart + readlen
            if startFeature < 0:
                startFeature = 0
            for (start, stop, gid, sense, ftype) in featList[startFeature:]:
                if tagStart > stop:
                    startFeature += 1
                    continue
                if start > stopPoint:
                    startFeature -= 100
                    break
                if sense == 'R':
                    sense = '-'
                else:
                    sense = '+'
                if start <= tagStart <= stop and sense == tSense:
                    try:
                        gidReadDict[gid].append(tagReadID)
                        if tagReadID in read2GidDict:
                            if gid not in read2GidDict[tagReadID]:
                                read2GidDict[tagReadID].append(gid)
                        else:
                            read2GidDict[tagReadID] = [gid]
                        gidCount[gid] += 1
                    except:
                        print "gid %s not in gidReadDict" % gid
                    stopPoint = stop

for gid in gidList:
    if 'FAR' not in gid:
        symbol = 'LOC' + gid
        geneinfo = ''
        try:
            geneinfo = geneinfoDict[gid]
            if genome == 'celegans':
                symbol = geneinfo[0][1]
            else:
                symbol = geneinfo[0][0]
        except:
            try:
                symbol = geneannotDict[(genome, gid)][0]
            except:
                symbol = 'LOC' + gid
                #print "problem with %s" % (gid)
    else:
        symbol = gid
    tagCount = 0.
    for readID in gidReadDict[gid]:
        try:
            tagValue = uniqueCountDict[gid]
        except:
            tagValue = 1
        tagDenom = 0.
        for aGid in read2GidDict[readID]:
            try:
                tagDenom += uniqueCountDict[aGid]
            except:
                tagDenom += 1
        tagCount += tagValue / tagDenom
    
    if doVerbose:
        print "%s %s %f" % (gid, symbol, tagCount)
    outfile.write("%s\t%s\t%d\n" % (gid, symbol, tagCount))
#infile.close()
outfile.close()


if '-cache' in sys.argv:
    uncacheGeneDB(genome)
