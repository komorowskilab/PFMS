try:
	import psyco
	psyco.full()
except:
	print 'psyco not running'

import sys

print '%s: version 5.1' % sys.argv[0]

if len(sys.argv) < 4:
	print 'usage: python %s genome rdsfile outfilename [-stranded] [-splices] [-noUniqs] [-multi] [-models file] [-replacemodels] [-searchGID] [-countfeatures] [-cache pages] [-markGID]' % sys.argv[0]
	sys.exit(1)

from commoncode import *
from cistematic.genomes import Genome
from cistematic.core.geneinfo import geneinfoDB

genome = sys.argv[1]
hitfile =  sys.argv[2]
outfilename = sys.argv[3]

searchGID = False
if '-searchGID' in sys.argv:
    searchGID = True

countFeats = False
if '-countfeatures' in sys.argv:
    countFeats = True

doStranded = 'both'
doUniqs = True
doMulti = False
doSplices = False
if '-stranded' in sys.argv:
    print "will track strandedness"
    doStranded = 'track'
if '-noUniqs' in sys.argv:
    doUniqs = False
if '-multi' in sys.argv:
    doMulti = True
if '-splices' in sys.argv:
    doSplices = True

extendGenome = ''
replaceModels = False
if '-models' in sys.argv:
    extendGenome = sys.argv[sys.argv.index('-models') + 1]
    if '-replacemodels' in sys.argv:
        replaceModels = True
        print "will replace gene models with %s" % extendGenome
    else:
        print "will extend gene models with %s" % extendGenome

markGID = False
if '-markGID' in sys.argv:
    print "marking GID and NM"
    markGID = True
    
doCache = False
cachePages = 100000
if '-cache' in sys.argv:
    doCache = True
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])

hitRDS = readDataset(hitfile, verbose = True, cache=doCache)
if cachePages > hitRDS.getDefaultCacheSize():
    hitRDS.setDBcache(cachePages)

readlen = hitRDS.getReadSize()

hg = Genome(genome, inRAM=True)
if extendGenome != '':
    hg.extendFeatures(extendGenome, replace = replaceModels)

idb = geneinfoDB(cache=True)

gidDict = {}
geneinfoDict = idb.getallGeneInfo(genome)
geneannotDict = hg.allAnnotInfo()
gidCount = {}

print "getting gene features...."
featuresByChromDict = getFeaturesByChromDict(hg)

seenFeaturesByChromDict = []
print "getting geneIDs...."
gidList = hg.allGIDs()

gidList.sort()
for gid in gidList:
    gidCount[gid] = 0

index = 0
chromList = hitRDS.getChromosomes(fullChrom=False)
if len(chromList) == 0 and doSplices:
    chromList = hitRDS.getChromosomes(table='splices', fullChrom=False)

if markGID:
    print "Flagging all reads as NM"
    hitRDS.setFlags('NM', uniqs=doUniqs, multi=doMulti, splices=doSplices)
    
for chrom in chromList:
    if chrom not in featuresByChromDict:
        continue
    if countFeats:
        seenFeaturesByChromDict[chrom] = []
    print '\nchr' + chrom
    startFeature = 0
    fullchrom = 'chr' + chrom
    regionList = []
    #if markGID:
    #    hitRDS.setSynchronousPragma('OFF')        
    print "counting GIDs"
    for (start, stop, gid, fsense, ftype) in featuresByChromDict[chrom]:
        try:
            if doStranded == 'track':
                checkSense = '+'
                if fsense == 'R':
                    checkSense = '-'
                regionList.append((gid, fullchrom, start, stop, checkSense))
                count = hitRDS.getCounts(fullchrom, start, stop, uniqs=doUniqs, multi=doMulti, splices=doSplices, sense=checkSense)
            else:
                regionList.append((gid, fullchrom, start, stop))
                count = hitRDS.getCounts(fullchrom, start, stop, uniqs=doUniqs, multi=doMulti, splices=doSplices)
            gidCount[gid] += count
            if countFeats:
                if (start, stop, gid, fsense) not in seenFeaturesByChromDict[chrom]:
                    seenFeaturesByChromDict[chrom].append((start, stop, gid, fsense))
        except:
            print "problem with %s - skipping" % gid
    if markGID:
        print "marking GIDs"
        hitRDS.flagReads(regionList, uniqs=doUniqs, multi=doMulti, splices=doSplices, sense=doStranded)
        #hitRDS.setSynchronousPragma('ON')
        print "finished marking"

#infile = open(infilename)
print ' '
if countFeats:
    count = 0
    for chrom in seenFeaturesByChromDict:
        count += len(seenFeaturesByChromDict[chrom])
    print "saw %d features" % count

outfile = open(outfilename,'w')

for gid in gidList:
    symbol = 'LOC' + gid
    geneinfo = ''
    thegid = gid
    if searchGID and gid not in geneinfoDict:
        actualGeneID = idb.getGeneID(genome, gid)
        if len(actualGeneID) > 0:
            thegid = actualGeneID[1]
    try:
        geneinfo = geneinfoDict[thegid]
        symbol = geneinfo[0][0]
        #print geneinfo
    except:
        try:
            symbol = geneannotDict[(genome, gid)][0]
        except:
            symbol = 'LOC' + gid
            #print "problem with %s" % (gid)
    if gid in gidCount:
        outfile.write("%s\t%s\t%d\n" % (gid, symbol, gidCount[gid]))
    else:
        outfile.write("%s\t%s\t0\n" % (gid, symbol))
if markGID and doCache:
    hitRDS.saveCacheDB(hitfile)
#infile.close()
outfile.close()

