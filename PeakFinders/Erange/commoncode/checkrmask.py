try:
    import psyco
    psyco.full()
except:
    pass

import sqlite3 as sqlite
import sys, string
import os.path
from commoncode import writeLog

versionString = '%s: version 3.5' % sys.argv[0]
print versionString

if len(sys.argv) < 4:
    print 'usage: python %s dbfile infile outfile goodfile [-startField field] [-cache numPages] [-log logfile]' % sys.argv[0]
    sys.exit(1)

dbfile = sys.argv[1]
filename = sys.argv[2]
startField = 0
outfile = open(sys.argv[3],'w')
goodfile = open(sys.argv[4],'w')
if '-startField' in sys.argv:
    startField = int(sys.argv[sys.argv.index('-startField') + 1])
if startField < 0:
    startField = 0

cachePages = 500000
if '-cache' in sys.argv:
    cachePages =  int(sys.argv[sys.argv.index('-cache') + 1])
    if cachePages < 250000:
        cachePages = 250000

doLog = False
if '-log' in sys.argv:
    logfilename = sys.argv[sys.argv.index('-log') + 1]
    writeLog(logfilename, versionString, string.join(sys.argv[1:]))
    doLog = True

infile = open(filename)
if os.path.isfile(dbfile):
    db = sqlite.connect(dbfile)
    sql = db.cursor()
    sql.execute("PRAGMA CACHE_SIZE = %d" % cachePages)
    sql.execute("PRAGMA temp_store = MEMORY")
else:
    print "No database - passing through"
    if doLog:
        writeLog(logfilename, versionString, "No database - passing through")    
    for line in infile:
        outfile.write(line + '\tNR\tNR\t0.00\n')
        goodfile.write(line)
    outfile.close()
    goodfile.close()
    sys.exit(0)

featureList = []
featureDict = {}

for line in infile:
    if line[0] == '#':
        continue
    fields = line.strip().split('\t')
    chrom = fields[startField][3:]
    start = int(fields[startField + 1])
    stop = int(fields[startField + 2])
    featureList.append((chrom,start, stop))
    featureDict[(chrom, start, stop)] = line.strip()
infile.close()

featureList.sort()
currentChrom=''
increment = 20000000
for (chrom, start, stop) in featureList:
	if chrom != currentChrom:
		currentMax = 0
	if start > currentMax:
		currentChrom = chrom
                currentMin = currentMax
                currentMax += increment
                print 'caching %s from %d to %d' % (chrom, currentMin, currentMax)
		try:
			del con
		except:
			pass
		con = sqlite.connect(":memory:")
		sql.execute('select start, stop, name, family from repeats where chrom = "%s" and start >= %d and start <= %d order by start' % (chrom, currentMin, currentMax + 10000))
		results = sql.fetchall()
		results2 = []
		con.execute("create table repeats(name, family, start, stop)")
                con.execute("PRAGMA CACHE_SIZE = %d" % cachePages)
                con.execute("PRAGMA temp_store = MEMORY")
		for (rstart, rstop, name, family) in results:
			results2.append((name, family, int(rstart), int(rstop)))
		con.executemany("insert into repeats(name, family, start, stop) values (?,?,?,?)", results2)
                con.execute("CREATE INDEX posIndex on repeats(start, stop)")
		print chrom, len(results2)
		sql2 = con.cursor()
	featureLength = abs(stop - start)
	results = []
	finalresults = []
	sql2.execute('select start, stop, name, family from repeats where start < %d and stop > %d' % (start, start))
	results = sql2.fetchall()
	for (rstart, rstop, name, family) in results:
		overlapLength = float(abs(rstop - start))
		if overlapLength > featureLength:
			overlapLength = featureLength
		ratio = overlapLength / featureLength
		if (name, family, ratio) not in finalresults:
			finalresults.append((name, family, ratio))
	sql2.execute('select start, stop, name, family from repeats where start < %d and stop > %d' % (stop, stop))
	results = sql2.fetchall()
	for (rstart, rstop, name, family) in results:
		overlapLength = float(abs(rstart - stop))
		if overlapLength > featureLength:
			overlapLength = featureLength
		ratio = overlapLength / featureLength
		if (name, family, ratio) not in finalresults:
			finalresults.append((name, family, ratio))
	sql2.execute('select start, stop, name, family from repeats where start <= %d and stop >= %d' % (start, stop)) 
	results = sql2.fetchall()
	for (rstart, rstop, name, family) in results:
		overlapLength = float(abs(rstop - rstart))
		if overlapLength > featureLength:
			overlapLength = featureLength
		ratio = overlapLength / featureLength
		if (name, family, ratio) not in finalresults:
			finalresults.append((name, family, ratio))
	sql2.execute('select start, stop, name, family from repeats where start >= %d and stop <= %d' % (start, stop)) 
	results = sql2.fetchall()
	for (rstart, rstop, name, family) in results:
		overlapLength = float(abs(rstop - rstart))
		if overlapLength > featureLength:
			overlapLength = featureLength
		ratio = overlapLength / featureLength
		if (name, family, ratio) not in finalresults:
			finalresults.append((name, family, ratio))
	line = featureDict[(chrom, start, stop)]
	total = 0.
	for (name, family, fraction) in finalresults:
		outline = '%s\t%s\t%s\t%2.2f' % (line, name, family, fraction)
		total += fraction
		print outline
		outfile.write(outline + '\n')
	if len(finalresults) == 0:
		outline = '%s\t%s\t%s\t%2.2f' % (line, 'NR', 'NR', 0.)
		print outline
		outfile.write(outline + '\n')
	if total < 0.2:
		goodfile.write(line + '\n')


