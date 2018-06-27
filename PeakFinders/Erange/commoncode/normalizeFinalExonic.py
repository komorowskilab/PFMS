try:
	import psyco
	psyco.full()
except:
	pass

import sys
print "%s: version 3.5" % sys.argv[0]
if len(sys.argv) < 5:
	print "usage: python %s rdsfile expandedRPKMfile multicountfile outfile [-multifraction] [-multifold] [-minrpkm minThreshold] [-cache] [-withGID]\n" % sys.argv[0]
	sys.exit(1)

from commoncode import *

rdsfilename = sys.argv[1]
expandedRPKMfile = open(sys.argv[2])
multicountfile = open(sys.argv[3])

writeGID = False
reportFraction = False
reportFold = False 
minThreshold = 0.

if '-withGID' in sys.argv:
    writeGID = True
    
if '-multifraction' in sys.argv:
    reportFraction = True
    print "reporting fractional contribution of multireads"
elif '-multifold' in sys.argv:
    reportFold = True
    print "reporting fold contribution of multireads"

if '-minrpkm' in sys.argv:
    minThreshold = float(sys.argv[sys.argv.index('-minrpkm') + 1])
    print "returning only genes with > %.2f RPKM" % minThreshold
    
doCache = False
if '-cache' in sys.argv:
    doCache = True

RDS = readDataset(rdsfilename, verbose = True, cache=doCache, reportCount=False)
uniqcount = RDS.getUniqsCount()
splicecount = RDS.getSplicesCount()
multicount = RDS.getMultiCount()
countDict = {}
multicountDict = {}
lengthDict = {}
gidList = []

uniqspliceCount = (uniqcount + splicecount) / 1000000.
totalCount = (uniqcount + splicecount + multicount) / 1000000.

symbolDict = {}

for line in expandedRPKMfile:
    fields = line.strip().split()
    lineGID = fields[0]
    symbolDict[lineGID] = fields[1]
    countDict[lineGID] = float(fields[-1]) * float(fields[-2]) * uniqspliceCount
    lengthDict[lineGID] = float(fields[-2])
    multicountDict[lineGID] = 0
    if lineGID not in gidList:
        gidList.append(lineGID)
expandedRPKMfile.close()

for line in multicountfile:
    fields = line.strip().split()
    gid = fields[0]
    if gid in countDict:
        countDict[gid] += float(fields[-1])
        multicountDict[gid] = float(fields[-1])
    else:
        print "could not find gid %s in dictionaries" % gid
        
multicountfile.close()

outfile = open(sys.argv[4],'w')
outheader = '#'
if writeGID:
    outheader += 'GID\t'
outheader += 'gene\tlen_kb\tRPKM'
if reportFraction:
    outheader += '\tmulti/all'
elif reportFold:
    outheader += '\tall/uniq'
outheader += '\n'
outfile.write(outheader)

outlineList = []
index = 0
for gid in gidList:
    outline = ''
    gene = symbolDict[gid]
    #rpkm is rpm / geneLength in kb.
    rpm = countDict[gid] / totalCount
    rpkm = rpm / lengthDict[gid]
    if rpkm < minThreshold:
        continue
    if writeGID:
        outline = '%s\t' % gid
    index += 1
    try:
        multirpm = multicountDict[gid] / totalCount
        multirpkm = multirpm / lengthDict[gid]
    except:
        print "problem with %s - skipping " % gid
        continue
    if reportFraction or reportFold:
        try:
            if reportFraction:
                multivalue = multirpkm / rpkm
            else:
                if rpm > multirpm:
                    uniqrpkm = (rpm - multirpm) / lengthDict[gid]
                    multivalue = rpkm / uniqrpkm
                elif rpkm > 0.01:
                    multivalue = 100.
                else:
                    multivalue = 1.0
        except:
            multivalue = 0
        outline += '%s\t%.3f\t%.2f\t%.2f\n' %  (gene, lengthDict[gid], rpkm, multivalue)
        outlineList.append((rpkm, outline))
    else:
        outline += '%s\t%.3f\t%.2f\n' %  (gene, lengthDict[gid], rpkm)
        outlineList.append((rpkm, outline))

outlineList.sort()
outlineList.reverse()

for (rpkm, line) in outlineList:
    outfile.write(line)
outfile.close()

print "returned %d genes" % index
