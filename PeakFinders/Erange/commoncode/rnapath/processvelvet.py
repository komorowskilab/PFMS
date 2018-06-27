import sys
print '%s: version 1.1' % sys.argv[0]

if len(sys.argv) < 3:
	print 'usage: python %s infile outfile [-prefix contigpref] [-filter pslfile] [-min bp] [-keepcov]' % sys.argv[0]
	sys.exit(1)

infile = open(sys.argv[1])
outfile = open(sys.argv[2],'w')

filterList = []
acceptedSize = 0
nSize = 0
contigsAccepted = 0
minSize = 0

keepCoverage = False
if '-keepcov' in sys.argv:
	keepCoverage = True
contigPrefix = 'chr'
if '-prefix' in sys.argv:
	contigPrefix = sys.argv[sys.argv.index('-prefix') + 1]
if '-min' in sys.argv:
	minSize = int(sys.argv[sys.argv.index('-min') + 1])

filteredSize = 0
if '-filter' in sys.argv:
	filterFile = open(sys.argv[sys.argv.index('-filter') + 1])
	for line in filterFile:
		if 'NODE' in line:
			fields = line.strip().split()
			exclude = fields[9]
			if exclude not in filterList:
				filterList.append(exclude)
	filterFile.close()

completeID = ''
currentSeq = ''
for line in infile:
	if '>NODE' in line:
		if len(completeID) > 5 and completeID not in filterList:
			fields = completeID.split('_')
			newID = fields[1]
			if keepCoverage:
				newID = fields[1] + '_' + fields[-1].strip()
			if len(currentSeq) >= minSize:
				outfile.write('>%s%s\n%s' % (contigPrefix, newID, currentSeq))
				acceptedSize += len(currentSeq) - currentSeq.count('\n')
				nSize += currentSeq.count('N')
				contigsAccepted += 1
		else:
			filteredSize += len(currentSeq) - currentSeq.count('\n')
		completeID = line.strip()[1:]
		currentSeq = ''
	else:
		currentSeq += line
if len(completeID) > 5 and completeID not in filterList:
	fields = completeID.split('_')
	newID = fields[1]
	if keepCoverage:
		newID = fields[1] + '_' + fields[-1].strip()
	outfile.write('>%s%s\n%s' % (contigPrefix, newID, currentSeq))
	acceptedSize += len(currentSeq) - currentSeq.count('\n')
	nSize += currentSeq.count('N')
	contigsAccepted += 1
infile.close()
outfile.close()

print '%d contigs accepted' % contigsAccepted
print '%d bp original' % (acceptedSize + filteredSize)
print '%d bp accepted' % acceptedSize
print '%d bp accepted N' % nSize
print '%d bp filtered\n' % filteredSize
