try:
	import psyco
	psyco.full()
except:
	pass
import sys
import sqlite3 as sqlite
import os

print 'version 2.0'
if len(sys.argv) < 3:
	print 'usage: python %s rmaskdir rmaskdbfile' % sys.argv[0]
	sys.exit(1)

rmaskdir = sys.argv[1]
rmaskdb = sys.argv[2]

files = os.listdir(rmaskdir)
db = sqlite.connect(rmaskdb)
sql = db.cursor()
sql.execute('create table repeats (chrom varchar, start int, stop int, name varchar, family varchar)')
sql.execute("PRAGMA temp_store = MEMORY")
sql.execute("PRAGMA DEFAULT_CACHE_SIZE = 500000")
db.commit()

for filename in files:
	if 'rmsk' not in filename:
		continue
	print filename
	infile = open(rmaskdir + '/' + filename)
	for entry in infile:
		fields = entry.strip().split('\t')
		chrom = fields[5][3:]
		start = int(fields[6])
		stop = int(fields[7])
		name = fields[10]
		family = fields[12]
		stmt = 'insert into repeats values("%s", %d, %d, "%s", "%s")' % (chrom, start, stop, name, family)
		sql.execute(stmt)
	db.commit()
print 'building index...'
sql.execute("PRAGMA SYNCHRONOUS = OFF")
sql.execute('create index chromIndex on repeats(chrom)')
sql.execute('create index mainIndex on repeats(chrom, start, stop)')
db.commit()
