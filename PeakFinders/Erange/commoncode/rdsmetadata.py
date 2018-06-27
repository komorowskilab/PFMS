#
#  rdsmetadata.py
#  ENRAGE
#
try:
    import psyco
    psyco.full()
except:
    pass

import sys
from commoncode import readDataset

print '%s: version 2.7' % sys.argv[0]
if len(sys.argv) < 2:
    print 'usage: python %s rdsfile [propertyName1::propertyValue1] ... [propertyNameN::propertyValueN] [-defaultcache size] [-index] [-dropindex] [-nocount] [-complexity] [-reset] [-initrna] [-cache pages]' % sys.argv[0]
    print '\nwhere the optional metadata name::value pairs are added to the existing dataset\n'
    sys.exit(1)

doCount = True
if '-nocount' in sys.argv:
    doCount = False

doCache = False
cachePages = -1
if '-cache' in sys.argv:
    doCache = True
    try:
        cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])
    except: 
        pass

datafile = sys.argv[1]
if '-initrna' in sys.argv:
    rds = readDataset(datafile, initialize=True, datasetType='RNA', verbose=True, cache=doCache)
else:
    rds = readDataset(datafile, verbose=True, reportCount=doCount, cache=doCache)

if cachePages > rds.getDefaultCacheSize():
    rds.setDBcache(cachePages)

cacheVal = 0
if '-defaultcache' in sys.argv:
    cacheVal = int(sys.argv[sys.argv.index('-defaultcache') + 1])
    rds.setDBcache(cacheVal, default=True)
    print "set default cache size to %d pages" % cacheVal

if '-reset' in sys.argv:
    print "clearing read flags"
    rds.resetFlags()

if '-dropindex' in sys.argv:
    try:
        rds.dropIndex()
    except:
        print 'could not drop index'

if '-index' in sys.argv:
    print "building index...."
    if cacheVal > 0:
        rds.buildIndex(cacheVal)
    else:
        rds.buildIndex()

if '-complexity' in sys.argv:
    print "calculating uniq read complexity..."
    uniqs = rds.getUniqsCount(distinct=False)
    distincts = rds.getUniqsCount(distinct=True)
    print '%d distincts / %d uniqs = %.2f' % (distincts, uniqs, float(distincts) / uniqs)

propertyList = []
for arg in sys.argv:
    if '::' in arg:
        (pname, pvalue) = arg.strip().split('::')
        print 'adding %s : %s' % (pname, pvalue)
        propertyList.append((pname, pvalue))
if len(propertyList) > 0:
    rds.insertMetadata(propertyList)



