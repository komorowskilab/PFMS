import sys, string

if len(sys.argv) < 4:
    print 'usage: python %s factorlabel bedinfilename regionoutfile' % sys.argv[0]
    sys.exit(1)

factor = sys.argv[1]
infilename = sys.argv[2]
outfilename = sys.argv[3]

index = 1
infile = open(infilename)
outfile = open(outfilename, 'w')
for line in infile:
    if 'track' in line:
        continue
    fields = line.split()
    line = string.join(fields, '\t')
    outfile.write('%s%d\t%s\n' % (factor, index, line))
    index += 1
infile.close()
outfile.close()

