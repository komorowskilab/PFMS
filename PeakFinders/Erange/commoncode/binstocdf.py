import sys

print 'version 1.0'
if len(sys.argv) < 3:
	print 'usage: python %s infile outfile' % sys.argv[0]
	sys.exit(1)

infilename = sys.argv[1]
outfilename = sys.argv[2]

infile = open(infilename)
outfile = open(outfilename, 'w')

for line in infile:
	fields = line.strip().split()
	if len(fields) < 4:
		continue
	total = int(fields[2])
	if total == 0:
		outfile.write(line)
		continue
	outfile.write('%s\t%s\t%s\t%s' % (fields[0], fields[1], fields[2], fields[3]))
	cum = 0
	for bin in fields[4:]:
		cum += int(bin)
		percent = 100 * cum / total
		outfile.write('\t%d' % percent)
	outfile.write('\n')
infile.close()
outfile.close()
