#
#  combinerds.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    pass

import sys
from commoncode import readDataset

print '%s: version 1.1' % sys.argv[0]
if len(sys.argv) < 2:
    print 'usage: python %s destinationRDS inputrds1 [inputrds2 ....] [-table table_name] [-init] [-initrna] [-index] [-cache pages]' % sys.argv[0]
    #print '\nwhere the optional metadata name::value pairs are added to the existing dataset\n'
    sys.exit(1)

doCache = False
cachePages = -1
if '-cache' in sys.argv:
    doCache = True
    try:
        cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])
    except: 
        pass

datafile = sys.argv[1]
infileList = []
for index in range(2, len(sys.argv)):
    if sys.argv[index][0] == '-':
        break
    infileList.append(sys.argv[index])

print "destination RDS: %s" % datafile

if '-initrna' in sys.argv:
    rds = readDataset(datafile, initialize=True, datasetType='RNA')
elif '-init' in sys.argv:
    rds = readDataset(datafile, initialize=True)

withFlag = ''
if '-flag' in sys.argv:
    withFlag = sys.argv[sys.argv.index('-flag') + 1]
    print "restrict to flag = %s" % withFlag

rds = readDataset(datafile, verbose=True, cache=doCache)

if cachePages > rds.getDefaultCacheSize():
    rds.setDBcache(cachePages)
    cacheVal = cachePages
else:
    cacheVal = rds.getDefaultCacheSize()

doIndex = False
if '-index' in sys.argv:
    doIndex = True

tableList = []
if '-table' in sys.argv:
    tableList.append(sys.argv[sys.argv.index('-table') + 1])
else:
    tableList = rds.getTables()

metaDict = rds.getMetadata()
if 'numberImports' not in metaDict:
    origIndex = 0
    rds.insertMetadata([('numberImports',str(0))])
else:
    origIndex = int(metaDict['numberImports'])

index = origIndex
for inputfile in infileList:
    asname = 'input' + str(index)
    rds.attachDB(inputfile,asname)
    for table in tableList:
        print "importing table %s from file %s" % (table, inputfile)
        ascols = '*'
        if table == 'uniqs':
            ascols = "NULL, '%s' || readID, chrom, start, stop, sense, weight, flag, mismatch" % asname
        elif table == 'multi':
            ascols = "NULL, '%s' || readID, chrom, start, stop, sense, weight, flag, mismatch" % asname
        elif table == 'splices':
            ascols = "NULL, '%s' || readID, chrom, startL, stopL, startR, stopR, sense, weight, flag, mismatch" % asname
        elif table == 'metadata':
            ascols = "name, value || ' (import_%d)'" % index
            rds.importFromDB(asname, table, ascols)
        
        if table != 'metadata':
            rds.importFromDB(asname, table, ascols, withFlag)
        
    rds.detachDB(asname)
    rds.insertMetadata([('import_' + str(index), '%s %s' % (inputfile, str(tableList)))])
    index += 1

rds.updateMetadata('numberImports', index, origIndex)
if doIndex:
    print "building index...."
    if cacheVal > 0:
        rds.buildIndex(cacheVal)
    else:
        rds.buildIndex()

if doCache:
    rds.saveCacheDB(datafile)
