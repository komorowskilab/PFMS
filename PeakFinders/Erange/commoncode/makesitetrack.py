#
#  makesitetrack.py
#  ENRAGE
#

import sys, string
print '%s: version 2.1' % sys.argv[0]

if len(sys.argv) < 3:
	print 'usage: python %s sitefile outbedfile [-noheader] [-stype fieldID] [-color xx,yy,zz] [-append] [-exploded]' % sys.argv[0]
	sys.exit(1)

color = '0,0,0'
stypeID = 4
doStype = False
if '-stype' in sys.argv:
    stypeID = int(sys.argv[sys.argv.index('-stype') + 1])
    doStype = True
    
if '-color' in sys.argv:
    color = sys.argv[sys.argv.index('-color') + 1]

infile = open(sys.argv[1])

if '-append' in sys.argv:
    outfile = open(sys.argv[2],'a')
else:
    outfile = open(sys.argv[2],'w')

compact = True
if '-exploded' in sys.argv:
    compact = False

try:
    (name, extension) = sys.argv[2].split('.')
except:
    name = sys.argv[2].split('.')[:-1]
    name = string.join(name, '_')
    
if '-noheader' not in sys.argv:
    outfile.write('track name="%s" visibility=4 itemRgb="On"\n' % name)

count = 1
for line in infile:
    if line[0] == '#':
        continue
    fields = line.split()
    if compact:
        (chrom, loc) = fields[0].split(':')
        (start, stop) = loc.split('-')
        score = fields[1]
    else:
        chrom = fields[1]
        start = fields[2]
        stop = fields[3]
        score = 1.
    stype = name + '-' + str(count)
    if doStype:
        try:
            stype = fields[stypeID]
            if stype == '11':
                stype = 'can'
            elif stype == '0':
                stype = 'half'
            else:
                stype = 'NC' + stype
        except:
            pass
    sense = fields[-2].strip()
    if sense not in ['+','-']:
        sense = '+'
    outfile.write('%s\t%s\t%d\t%s\t%s\t%s\t-\t-\t%s\n' % (chrom, start, int(stop) + 1, stype, score, sense, color))
    count += 1
outfile.close()
