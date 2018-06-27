#
#  peakstoregion.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
	pass
import sys

print '%s: version 1.0' % sys.argv[0]
if len(sys.argv) < 3:
    print 'usage: python %s peakfile outfile [radius] [chromField] [posField] [labelField] [datafield]' % sys.argv[0]
    sys.exit(1)

peakfile = open(sys.argv[1])
outfile = open(sys.argv[2],'w')

radius = 500
chromField = 2
posField = 3
labelField = 1
dataField = -1

if len(sys.argv) > 3:
    radius = int(sys.argv[3])

if len(sys.argv) > 4:
    chromField = int(sys.argv[4])

if len(sys.argv) > 5:
    posField = int(sys.argv[5])

if len(sys.argv) > 6:
    labelfield = int(sys.argv[6])

if len(sys.argv) > 7:
    dataField = int(sys.argv[7])
    
for line in peakfile:
    fields = line.strip().split()
    label = 'REGION'
    try:
        label = fields[labelField]
    except:
        pass
    start = int(fields[posField]) - radius
    stop = int(fields[posField]) + radius
    outfile.write('%s\t%s\t%d\t%d\t%s\n' % (label, fields[chromField], start, stop, fields[dataField]))
outfile.close()
