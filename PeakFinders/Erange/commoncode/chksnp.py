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

print 'version 3.6: %s' % sys.argv[0]
if len(sys.argv) < 4:
    print 'usage: python %s dbfile snpsfile dbsnp_outfile [-cache numPages]' % sys.argv[0]
    sys.exit(1)

dbfile = sys.argv[1]
filename = sys.argv[2]
outfile = sys.argv[3]

print 'using %s ' % dbfile
cachePages = 500000
doCache = False
if '-cache' in sys.argv:
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])
    if cachePages < 500000:
        cachePages = 500000
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

infile = open(filename)
featureList = []
featureDict = {}

for line in infile:
    if line[0] == '#':
        continue
    fields = line.strip().split('\t')
    chrom = fields[2][3:]
    pos = int(fields[3])
    #print chrom, pos
    featureList.append((chrom,pos))
    featureDict[(chrom, pos)] = line.strip()
featureList.sort()

index = 0
for (chrom, pos) in featureList:
    index += 1
    results = []
    sql.execute('select func, name from snp where chrom = "%s" and start = %d and stop = %d' % (chrom, pos-1, pos)) 
    results = sql.fetchall()
    isSnp = "\tN\A\tN\A"
    found = False
    for (func, name) in results:
        isSnp = "\t" + str(name) + "\t" + str(func)
        found =True
    if not found:
        sql.execute('select func, name from snp where chrom = "%s" and start <= %d and stop >= %d' % (chrom, pos-1, pos)) 
        results = sql.fetchall()
        isSnp = "\tN\A\tN\A"
        for (func, name) in results:
            isSnp = "\t" + str(name) + "\t" + str(func)
        
    featureDict[(chrom,pos)] += isSnp
    if index % 100 == 0:
        print '.',
        sys.stdout.flush()
        
if doCache:
    print '\nremoving cache'
    del db
    os.remove(cachefile)
outF = open(outfile, 'w') 
#outStr = "Sl" + "\t" + "Cl" + "\t" + "chrom" + "\t" + "mis pos" + "\t\t" + "match" + "\t" + "uniq_mis"+ "\t" + "tot_mis" + "\t" + "base_chg" + "\t" + "known_snp" + "\t" + "function" + "\n"
outStr = ''
outF.write(outStr)
for key,value in featureDict.iteritems():
    outStr = str(value) + "\n"
    outF.write(outStr)
outF.close()




