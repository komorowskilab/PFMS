#
#  buildMatrix.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 3/6/09.
#
import sys, string
from commoncode import writeLog

versionString = '%s: version 1.3' % sys.argv[0]
print versionString

if len(sys.argv) < 4:
    print 'usage: python %s matrix.step.N-1 data.part matrix.step.N [-rescale] [-truncate maxRPKM] [-log altlogfile]' % sys.argv[0]
    sys.exit(0)

infile = open(sys.argv[1])
colfilename = sys.argv[2]
outfilename = sys.argv[3]

rescale = False

logfilename = 'buildMatrix.log'
if '-log' in sys.argv:
    logfilename = sys.argv[sys.argv.index('-log') + 1]

truncateRPKM = False
maxRPKM = 100000000
if '-truncate' in sys.argv:
    truncateRPKM = True
    maxRPKM = int(sys.argv[sys.argv.index('-truncate') + 1])

if '-rescale' in sys.argv:
    rescale = True

writeLog(logfilename, versionString, string.join(sys.argv[1:]))

if '/' in colfilename:
    colname = colfilename.split('/')[-1]
else:
    colname = colfilename
fileParts = colname.split('.')
colID =  fileParts[0]

colfile = open(colfilename)
outfile = open(outfilename,'w')
header = infile.readline()[:-1]
if header.strip() == '':
    header = '#\t'

outfile.write(header + '\t' + colID +'\n')

values = []
min = 20000000000.
max = -1.
untruncatedMax = -1.
for line in colfile:
    if line[0] == '#':
        continue	
    fields = line.strip().split()
    val = float(fields[-1])
    if truncateRPKM and val > maxRPKM:
        if val > untruncatedMax:
            untruncatedMax = val
        val = maxRPKM
    values.append(val)
    if val < min:
        min = val
    if val > max:
        max = val

range = max - min
if rescale:
    finalValues = [(val - min)/range for val in values]
else:
    finalValues = values

for val in finalValues:
    line = infile.readline().strip()
    line += '\t%1.3f\n' % val
    outfile.write(line)
outfile.close()

if untruncatedMax > 0:
    max = untruncatedMax
message = "max value in %s was %.2f" % (colname, max)
if untruncatedMax > 0:
    message += " but was truncated to %d" % maxRPKM
print message
writeLog(logfilename, versionString, message)
