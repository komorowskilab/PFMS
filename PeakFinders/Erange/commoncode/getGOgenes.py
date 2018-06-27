import sys

from cistematic.genomes import Genome
from cistematic.core.geneinfo import geneinfoDB

print '%s: version 3.1' % sys.argv[0]

if len(sys.argv) < 3:
    print 'usage: python %s genome GOID1 [GOID2 ....] [-outfile outfilename] [-append] [-restrict genefile]'
    sys.exit(1)

genome = sys.argv[1]

writeOut = False
if '-outfile' in sys.argv:
    writeOut = True
    outfilename = sys.argv[sys.argv.index('-outfile') + 1]

restrict = False
if '-restrict' in sys.argv:
    restrictfilename = sys.argv[sys.argv.index('-restrict') + 1]
    restrict = True
    
hg = Genome(genome)
idb = geneinfoDB()

GOIDlist = []
for arg in sys.argv:
    if 'GO:' in arg:
        GOIDlist.append(arg)

print sys.argv
print GOIDlist

firstGeneList = []
for GOID in GOIDlist:
    testList = hg.allGIDsbyGOID(GOID)
    print 'GOID: %s (%d)' % (GOID, len(testList))
    firstGeneList += testList

geneDict = {}
for gid in firstGeneList:
    geneDict[gid] = 1

geneList = geneDict.keys()
    
print len(geneList)
geneInfoList = idb.getallGeneInfo(genome)

if writeOut:
    if '-append' in sys.argv:
        outfile = open(outfilename,'a')
    else:
        outfile = open(outfilename,'w')
    for GOID in GOIDlist:
        outfile.write('#%s\n' % GOID) 

restrictList = []
restrictDict = {}
symbolDict = []
if restrict:
    restrictFile = open(restrictfilename)
    for line in restrictFile:
        fields = line.strip().split()
        restrictList.append(fields[0])
        restrictDict[fields[0]] = line

outList = []
symbolDict = {}
for gid in geneList:
    symbol = 'LOC' + gid
    if restrict and gid not in restrictList:
        continue
    try:
        symbol = geneInfoList[gid][0][0]
    except:
        pass
    if restrict:
        symbolDict[symbol] = restrictDict[gid]
    outList.append(symbol)

outList.sort()
for symbol in outList:
    if writeOut:
        if restrict:
            outfile.write(symbolDict[symbol])
        else:
            outfile.write(symbol + '\n')
    else:
        print symbol

if writeOut:
    outfile.close()
