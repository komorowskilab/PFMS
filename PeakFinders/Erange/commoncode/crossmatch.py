try:
	import psyco
	psyco.full()
except:
	pass

import sys
from cistematic.core.orthomatcher import orthoMatcher

print 'version 1.1'
if len(sys.argv) < 7:
    print 'usage: python %s prefix directory genome1 genefile1 genome2 genefile2 [genome3 genefile3 .....]' % sys.argv[0]
    sys.exit(1)

prefix = sys.argv[1]
directory = sys.argv[2]
matchFiles = {}

genomesToMatch = (len(sys.argv) - 3) / 2
for index in range(genomesToMatch):
    genome = sys.argv[3 + index * 2]
    print genome
    if genome not in matchFiles:
        matchFiles[genome] = []
    matchFiles[genome].append(sys.argv[4 + index * 2])

print matchFiles
orthoMatcher(matchFiles, prefix, directory, fileList=True)
