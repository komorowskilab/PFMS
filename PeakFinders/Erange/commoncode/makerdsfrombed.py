#
#  makerdsfrombed.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 6/21/08.
#
try:
    import psyco
    psyco.full()
except:
    pass

import sys, string
from commoncode import readDataset, writeLog

verstring = "%s: version 2.1" % sys.argv[0]
print verstring

if len(sys.argv) < 3:
    print "usage: python %s label bedfile outrdsfile [-append] [-index] [propertyName::propertyValue] [-cache numPages]" % sys.argv[0]
    print "\ntreats all imported reads as uniquely mapped\n"
    sys.exit(1)

label = sys.argv[1]
filename = sys.argv[2]
outdbname = sys.argv[3]

init = True
if '-append' in sys.argv:
    init = False

dataType = 'DNA'
if '-RNA' in sys.argv:
    dataType = 'RNA'

doIndex = False
if '-index' in sys.argv:
    doIndex = True

cachePages = 100000
if '-cache' in sys.argv:
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])

readsize = 0
padsize = 0
index = 0
insertSize = 100000

writeLog(outdbname + '.log', verstring, string.join(sys.argv[1:]))

infile = open(filename,'r')
#infile.readline()

rds = readDataset(outdbname, init, dataType, verbose=True)
if not init:
    rds.dropIndex()

#check that our cacheSize is better than the dataset's default cache size
defaultCacheSize = rds.getDefaultCacheSize()
if cachePages > defaultCacheSize:
    if init:
        rds.setDBcache(cachePages, default=True)
    else:
        rds.setDBcache(cachePages)

propertyList = []
for arg in sys.argv:
    if '::' in arg:
        (pname, pvalue) = arg.strip().split('::')
        propertyList.append((pname, pvalue))
if len(propertyList) > 0:
    rds.insertMetadata(propertyList)

insertList = []
for line in infile:
    if 'track' in line:
        continue
    fields = line.split()
    if readsize == 0:
        readsize = abs(int(fields[1]) - int(fields[2]))
        if init:
            rds.insertMetadata([('readsize', readsize+1)])
            rds.insertMetadata([('imported_from_bed', 'True')])
    chrom = fields[0]
    start = int(fields[1])
    stop = int(fields[2])
    sense = fields[5]
    insertList.append((label + '-' + str(index), chrom, start, stop, sense, 1.0, '', ''))
    if index % insertSize == 0:
        rds.insertUniqs(insertList)
        insertList = []
        print '.',
        sys.stdout.flush()
    index += 1

if len(insertList) > 0:
    rds.insertUniqs(insertList)

countString = '%d unique reads' % index
print countString

writeLog(outdbname + '.log', verstring, countString)

if doIndex:
    print "building index...."
    if cachePages > defaultCacheSize:
        rds.setDBcache(cachePages)
        rds.buildIndex(cachePages)
    else:
        rds.buildIndex(defaultCacheSize)

