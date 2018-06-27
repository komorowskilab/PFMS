import sys

print 'version 1.2'
if len(sys.argv) < 3:
	print 'usage: python %s field filename' % sys.argv[0]
        print '\n\tlike all good python programs, fields are counted starting at zero.\n'
	sys.exit(1)

fieldID = int(sys.argv[1])
filename = sys.argv[2]

infile = open(filename)
count = 0

for line in infile:
	fields = line.strip().split()
	try:
		if '.' in fields[fieldID]:
			count += float(fields[fieldID])
		else:
			count += int(fields[fieldID])
	except:
		pass
infile.close()
print count
