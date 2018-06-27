try:
	import psyco
	psyco.full()
except:
	pass
import sys, math

print '%s: version 3.1' % sys.argv[0]
if len(sys.argv) < 4:
	print 'usage: python %s label regionfile outbedfile [-color r,g,b] [-score field] [-narrowPeak] [-broadPeak] [-itemRgb] [-nolabel]\n' % sys.argv[0]
	print '\twhere color is in comma-delimited RGB without space'
        print '\tand field is a column with a score (first column is 0, second is 1,...)\n'  
        print '\t-narrowPeak assumes that findall.py was run with -listPeak'
        print '\t-broadPeak assumes that findall.py was *NOT* run with -listPeak\n'
	sys.exit(1)

factorlabel = sys.argv[1]
regionfile = open(sys.argv[2])

color = '0,0,0'
if '-color' in sys.argv:
    color = sys.argv[sys.argv.index('-color') + 1]

scoreField = None
if '-score' in sys.argv:
    scoreField = int(sys.argv[sys.argv.index('-score') + 1])

doNarrow = False
if '-narrowPeak' in sys.argv:
    doNarrow = True

doBroad = False
if '-broadPeak' in sys.argv:
    doBroad = True

itemRGB = False
if '-itemRGB' in sys.argv:
    itemRGB = True
    print "assigning each item its color"

outfile = open(sys.argv[3], 'w')

if '-nolabel' not in sys.argv:
    if itemRGB:
        outfile.write('track name=%s visibility=4 itemRgb="on"\n' % factorlabel)
    else:
        outfile.write('track name=%s visibility=4 color=%s\n' % (factorlabel, color))

for line in regionfile:
    if line[0] == '#':
        continue
    fields = line.strip().split()
    if doNarrow:
        signalVal = float(fields[4])
        pval = float(fields[-1])
        if pval == 0.:
            pValue = 350
        else:
            pValue = -1. * math.log(pval, 10)
        peakPos = int(fields[9]) - int(fields[2])
        outfile.write('%s\t%s\t%s\t%s\t%d\t.\t%.4f\t%.4f\t-1\t%d' % (fields[1], fields[2], fields[3], fields[0], 0, signalVal, pValue, peakPos))
    elif doBroad:
        signalVal = float(fields[4])
        pval = float(fields[-1])
        if pval == 0.:
            pValue = 350
        else:
            pValue = -1. * math.log(pval, 10)
        outfile.write('%s\t%s\t%s\t%s\t%d\t.\t%.4f\t%.4f\t-1' % (fields[1], fields[2], fields[3], fields[0], 0, signalVal, pValue))
    elif scoreField:
        score = int(float(fields[scoreField]))
        if score > 1000:
            score = 1000
        outfile.write('%s\t%s\t%s\t%s\t%s' % (fields[1], fields[2], fields[3], fields[0], score))
        if itemRGB:
            outfile.write('\t+\t-\t-\t%s' % color)
    else:
        outfile.write('%s\t%s\t%s\t%s' % (fields[1], fields[2], fields[3], fields[0]))
        if itemRGB:
            outfile.write('\t1000\t+\t-\t-\t%s' % color)
    outfile.write('\n')
outfile.close()
