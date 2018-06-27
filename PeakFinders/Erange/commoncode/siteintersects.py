#
#  siteintersects.py
#  ENRAGE
#

import sys

print '%s: version 2.0' % sys.argv[0]
if len(sys.argv) < 4:
    print 'usage: python %s sitefile1 sitefile2 outfile [-reject rejectfile1 rejectfile2] [-expanded]' % sys.argv[0]
    sys.exit(1)

sitefilename1 =  sys.argv[1]
sitefilename2 = sys.argv[2]

outfilename = sys.argv[3]
outfile = open(outfilename, 'w')

doReject = False
if '-reject' in sys.argv:    
    reject1file = open(sys.argv[sys.argv.index('-reject') + 1],'w')
    reject2file = open(sys.argv[sys.argv.index('-reject') + 2],'w')
    doReject = True

doExpanded = False
if '-expanded' in sys.argv:
    doExpanded = True
siteDict = {}
file1Dict = {}
unique2List = []

infile1count = 0
infile = open(sitefilename1)
infile.readline()
for line in infile:
    if line[0] == '#':
        continue
    infile1count += 1
    fields = line.strip().split()
    if doExpanded:
        chrom = fields[1][3:]
        start = int(fields[2])
        stop = int(fields[3])
        rest = fields[4:]
    else:
        (chrom, pos) = fields[0].split(':')
        chrom = chrom[3:]
        (start, stop) = pos.split('-')
        start = int(start)
        stop = int(stop)
        rest = fields[1:]
    try:
        siteDict[chrom].append((start, stop, rest))
    except:
        siteDict[chrom] = [(start, stop, rest)]
    if doReject:
        file1Dict[str((chrom, start, stop, rest))] = line
infile.close()

print "file1: %d" % infile1count

infile2count = 0
infile = open(sitefilename2)
infile.readline()

commonSites = 0
for line in infile:
    if line[0] == '#':
        continue
    infile2count += 1
    fields = line.strip().split()
    if doExpanded:
        chrom = fields[1][3:]
        start = int(fields[2])
        stop = int(fields[3])
        rest = fields[4:]
    else:
        (chrom, pos) = fields[0].split(':')
        chrom = chrom[3:]
        (start, stop) = pos.split('-')
        rest = str(fields[1:])
    start = int(start)
    stop = int(stop)
    mid = start + abs(stop - start)/2
    if chrom not in siteDict:
        if doReject:
            unique2List.append(line)
            continue
    twoNotCommon = True
    for (rstart, rstop, rline) in siteDict[chrom]:
        rsize = abs(rstart - rstop) /2
        rmid = rstart + abs(rstop - rstart)/2
        if abs(mid - rmid) < rsize:
            commonSites += 1
            if twoNotCommon:
                outfile.write('common%d\tchr%s\t%d\t%d\t%s\tchr%s\t%d\t%d\t%s\n' % (commonSites, chrom, rstart, rstop, str(rline), chrom, start, stop, rest))
                twoNotCommon = False
            try:
                if doReject:
                    del file1Dict[str((chrom, rstart, rstop, rline))]
            except:
                pass
                #print 'problem deleting %s' % str((chrom, rstart, rstop, rline))
    if doReject and twoNotCommon:
        unique2List.append(line)
outfile.close()

print "file2: %d" % infile2count

if doReject:
    for key in file1Dict:
        reject1file.write(file1Dict[key])
    for line in unique2List:
        reject2file.write(line)
    reject1file.close()
    reject2file.close()
print commonSites

