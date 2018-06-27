import sys

if len(sys.argv) < 4:
	print 'usage: python %s bins percent infile' % sys.argv[0]
	sys.exit(1)

bins = int(sys.argv[1])
percent = int(sys.argv[2])
infilename = sys.argv[3]

infile = open(infilename)
binsList = [0] * bins

for line in infile:
	fields = line.strip().split()
	index = 0
	for binCdf in fields[-1 * bins:]:
		if int(binCdf) > percent:
			binsList[index] += 1
			break
		index += 1
infile.close()
print binsList
