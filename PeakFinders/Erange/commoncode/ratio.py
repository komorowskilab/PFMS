import sys
import math
if len(sys.argv) < 3:
    print '%s: version 2.3' % sys.argv[0]
    print 'usage: python %s denominatorField infile [-only fieldID] [-out outfile]' % sys.argv[0]
    sys.exit(1)

field = int(sys.argv[1])
if sys.argv[2] == '-':
    infile = sys.stdin
else:
    infile = open(sys.argv[2])
record = False
if '-out' in sys.argv:
    print '%s: version 2.3' % sys.argv[0]
    outfile = open(sys.argv[sys.argv.index('-out') + 1],'w')
    record = True

doOnly = False
onlyField = -1
if '-only' in sys.argv:
    doOnly = True
    onlyField = int(sys.argv[sys.argv.index('-only') + 1])

line = infile.readline()
count = len(line.strip().split())
if record:
    outfile.write(line)
for line in infile:
    fields = line.strip().split()
    outline = '%s ' % fields[0]
    outError = False
    for index in range(1, count):
        if field == index:
            outline += '0 '
        elif doOnly and index != onlyField:
            outline += '%s ' % fields[index]
        else:
            try:
                outline += '%2.2f ' % math.log((float(fields[index]) + 1)/(float(fields[field]) + 1), 2)
            except:
                try:
                    outline += 'e%s ' % fields[index]
                except:
                    outError = True
    if outError:
        continue
    if record:
        outfile.write(outline + '\n')
    else:
        print outline
