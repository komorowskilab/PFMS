#
#  gointersects.py
#  ENRAGE
#

import sys

print '%s: version 1.0' % sys.argv[0]

if len(sys.argv) < 4:
    print 'usage: python %s gogidfile gidfile outfile' % sys.argv[0]
    sys.exit(1)

gogidfilename = sys.argv[1]
gidfilename = sys.argv[2]
outfilename = sys.argv[3]

gidList = []
gogidfile = open(gogidfilename)
for line in gogidfile:
	fields = line.split()
	gidList.append(fields[0])
gogidfile.close()

gidfile = open(gidfilename)
outfile = open(outfilename, 'w')
for line in gidfile:
	fields = line.split()
	if fields[0] in gidList:
		outfile.write(line)
gidfile.close()
outfile.close()
