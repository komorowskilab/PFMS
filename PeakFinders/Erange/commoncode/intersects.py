#
#  intersects.py
#  ENRAGE
#

import sys

print 'version 2.0'
if len(sys.argv) < 4:
    print 'usage: python %s infile1 infile2 [-d delimiter] [-file3 infile3] [-1 matchfield2] [-2 matchfield2] [-3 matchfield3] [-reject1 rejectfile] [-trackGID] outfile\n' % sys.argv[0]
    sys.exit(1)

infile1 = open(sys.argv[1])
infile2 = open(sys.argv[2])

matchField1 = 0
matchField2 = 0

delimiter = '\t'
doFile3 = False
trackGID = False
matchField3 = 0
if '-file3' in sys.argv:
    infile3 = open(sys.argv[sys.argv.index('-file3') + 1])
    doFile3 = True
    if '-3' in sys.argv:
        matchField3 = int(sys.argv[sys.argv.index('-3') + 1])

if '-1' in sys.argv:
    matchField1 = int(sys.argv[sys.argv.index('-1') + 1])

if '-2' in sys.argv:
    matchField2 = int(sys.argv[sys.argv.index('-2') + 1])

if '-d' in sys.argv:
    delimiter = sys.argv[sys.argv.index('-d') + 1]

doReject1 = False
if '-reject1' in sys.argv:
    doReject1 = True
    reject1file = open(sys.argv[sys.argv.index('-reject1') + 1],'w')

gidDict = {}
if '-trackGID' in sys.argv:
    trackGID = True

list1 = []
list2 = []
list3 = []
matchedList = []
matchedList12 = []
matchedList13 = []
matchedList23 = []

for line in infile1:
    if line[0] == '#':
        continue
    fields = line.strip().split(delimiter)
    candidate = fields[matchField1]
    if candidate not in list1:
        list1.append(candidate)
    if trackGID and candidate not in gidDict:
        gidDict[candidate] = fields[matchField1 + 1]

for line in infile2:
    if line[0] == '#':
        continue
    fields = line.strip().split(delimiter)
    candidate = fields[matchField2]
    if candidate not in list2:
        list2.append(candidate)
    if trackGID and candidate not in gidDict:
        gidDict[candidate] = fields[matchField2 + 1]

if doFile3:
    for line in infile3:
        if line[0] == '#':
            continue
        fields = line.strip().split(delimiter)
        candidate = fields[matchField3]
        if candidate not in list3:
            list3.append(candidate)
        if trackGID and candidate not in gidDict:
            gidDict[candidate] = fields[matchField3 + 1]
    infile3.close()

infile1.close()
infile2.close()

for candidate in list1:
    if doFile3 and candidate in list2 and candidate in list3:
        matchedList.append(candidate)
    elif doFile3 and candidate in list3:
        matchedList13.append(candidate)
    elif doFile3 and candidate in list2:
        matchedList12.append(candidate)
    elif not doFile3 and candidate in list2:
        matchedList.append(candidate)
    elif doReject1:
        if trackGID:
            reject1file.write('%s%s%s\n' % (candidate, delimiter, gidDict[candidate]))
        else:
            reject1file.write('%s\n' % candidate)

if doFile3:
    for candidate in list2:
        if candidate not in list1 and candidate in list3:
            matchedList23.append(candidate)

print len(list1), len(list2), len(list3)
if doFile3:
    print len(matchedList12), len(matchedList13), len(matchedList23)
print len(matchedList)

outfile = open(sys.argv[-1],'w')
for match in matchedList:
    if trackGID:
        outfile.write('%s%s%s\n' % (match, delimiter, gidDict[match]))
    else:
        outfile.write('%s\n' % match)

outfile.close()