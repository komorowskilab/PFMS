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
	print 'usage: python %s snpfile snpdbname' % sys.argv[0]
	sys.exit(1)

snpfilename = sys.argv[1]
snpdb = sys.argv[2]

db = sqlite.connect(snpdb)
sql = db.cursor()
sql.execute('create table snp (chrom varchar, start long, stop long, name varchar, observed varchar, strand varchar, ucscref varchar, ncbiref varchar, func varchar, moltype varchar, valid varchar, class varchar)')
sql.execute("PRAGMA temp_store = MEMORY")
sql.execute("PRAGMA DEFAULT_CACHE_SIZE = 500000")
db.commit()
#create table snp (chrom varchar, start int, stop int, name varchar, observed varchar, strand varchar, ucscref varchar, ncbiref varchar, func varchar, moltype varchar, valid varchar, class varchar)
# sample line in dbsnp file
# 608   chr1    3093453 3093454 rs52602943      0       +       G       G        C/G   genomic single  unknown 0       0       unknown exact   1

insertSize = 100000
insertCounter = 0
valuesList = []
print snpfilename
infile = open(snpfilename)
for entry in infile:
    try:
        fields = entry.strip().split('\t')
        chrom = fields[1][3:]
        start = int(fields[2])
        stop = int(fields[3])
        name = fields[4]
        strand = fields[6]
        refNcbi = fields[7]
        refUcsc = fields[8]
        observed = fields[9]
        molType = fields[10]
        classes = fields[11]
        valid = fields[12]
        func = fields[15]
    
        valuesList.append((chrom, start, stop, name, observed, strand, refUcsc, refNcbi, func, molType, valid, classes))
        insertCounter += 1
    except:
        continue
    if insertCounter % insertSize == 0:
        print insertCounter
        db.executemany("insert into snp values (?,?,?,?,?,?,?,?,?,?,?,?)", valuesList)
        valuesList = []
if len(valuesList) > 0:
    db.executemany("insert into snp values (?,?,?,?,?,?,?,?,?,?,?,?)", valuesList)
db.commit()

print 'building index'
sql.execute("PRAGMA SYNCHRONOUS = OFF")
sql.execute('create index chromIndex on snp(chrom)')
sql.execute('create index mainIndex on snp(chrom,start,stop)')
db.commit()
