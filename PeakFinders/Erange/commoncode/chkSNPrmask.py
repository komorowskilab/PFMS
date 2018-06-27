try:
    import psyco
    psyco.full()
except:
    pass

import sqlite3 as sqlite
import sys

import tempfile, shutil, os
from os import environ

if environ.get("CISTEMATIC_TEMP"):
    cisTemp = environ.get("CISTEMATIC_TEMP")
else:
    cisTemp = '/tmp'
tempfile.tempdir = cisTemp

print 'version 3.3: %s' % sys.argv[0]
if len(sys.argv) < 4:
    print 'usage: python %s dbfile snpsfile nr_snps_outfile [-cache numPages] [-repeats]' % sys.argv[0]
    sys.exit(1)

dbfile = sys.argv[1]
filename = sys.argv[2]
outfile = sys.argv[3]
repeats = False

print dbfile

if '-repeats' in sys.argv:
    repeats = True

cachePages = 500000
doCache = False
if '-cache' in sys.argv:
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])
    if cachePages < 250000:
        cachePages = 250000
    print 'caching locally...'
    cachefile = tempfile.mktemp() + '.db'
    shutil.copyfile(dbfile, cachefile)
    db = sqlite.connect(cachefile)
    doCache = True
    print 'cached...'
else:
    db = sqlite.connect(dbfile)

sql = db.cursor()
sql.execute("PRAGMA CACHE_SIZE = %d" % cachePages)
sql.execute("PRAGMA temp_store = MEMORY")
sql.execute("ANALYZE")

infile = open(filename)
featureList = []
featureDict = {}

for line in infile:
    if line[0] == '#':
        continue
    fields = line.strip().split('\t')
    chrom = fields[2][3:]
    pos = int(fields[3])
    featureList.append((chrom,pos))
    featureDict[(chrom, pos)] = line.strip()
featureList.sort()

index = 0
currentChrom=None
for (chrom, pos) in featureList:
    index += 1
    if chrom != currentChrom:
        print '\n%s' % chrom
        currentChrom = chrom
    results = []
    try:
        sql.execute('select family from repeats where chrom = %s and %d between start and stop' % (chrom, pos)) 
        results = sql.fetchall()
    except:
        pass
    if repeats: # if user wants to keep track of the SNPs in repeats
        featureDict[(chrom,pos)] += "\tN\A" 
        for x in results:
            featureDict[(chrom,pos)] += "\t" + str(x)
    else:
        for x in results:
            try:
                del featureDict[(chrom,pos)]
            except:
                pass
    if index % 100 == 0:
        print '.',
        sys.stdout.flush()

if doCache:
    print 'removing cache'
    del db
    os.remove(cachefile)
outF = open(outfile, 'w') 
for key,value in featureDict.iteritems():
    outStr = str(value) + "\n"
    outF.write(outStr)
outF.close()




