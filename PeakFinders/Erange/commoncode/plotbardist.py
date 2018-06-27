#
#  plotbardist.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 12/13/07.

import sys
try:
    import psyco
    psyco.full()
except:
	pass

def errorExit():
    print "usage: python %s infile1 [infile2] [infile3] [-bins numBins] [-field fieldNum] [-binSize num] [-doLog base] [-ymax maxY] [-xlabel label] [-ylabel label] [-binLabels labelList] [-title figtitle] [-legend legendList] [-xoffset value] [-figsize x,y] outfile.png" % sys.argv[0]
    print "where labelList and legendList are comma delimited strings of the form 'labelA,labelB,...,labelN'"
    sys.exit(1)

print '%s: version 3.2' % (sys.argv[0])
if len(sys.argv) < 3:
    errorExit()

numArgs = 0
for arg in sys.argv:
    if arg[0] == '-':
        # check that this isn't a negative number
        try:
            float(arg)
        except:
            numArgs += 1
numFiles = len(sys.argv) - 2 - 2 * numArgs

if numFiles < 1 or numFiles > 3:
    errorExit()

import matplotlib
matplotlib.use('Agg')

from pylab import *
from math import *

infile = open(sys.argv[1])
pngfilename = sys.argv[-1]

bins = 10
binnedField = -1
binLength = -1
doLog = False
logBase = 10
maxY = 0
yLabel = 'count'
xLabel = 'bins'
figTitle = ''
width = 0.5
offset = [-0.25]
colorList = ['b', 'r', 'c']

fileList = [infile]
if numFiles > 1:
    infile2 = open(sys.argv[2])
    fileList = [infile, infile2]
    width = 0.3
    offset = [-0.3,0]
if numFiles > 2:
    infile3 = open(sys.argv[3])
    fileList = [infile, infile2, infile3]
    width = 0.2
    offset = [-0.2,0.,0.2]
    
if '-xlabel' in sys.argv:
    xLabelField = sys.argv.index('-xlabel') + 1
    xLabel = sys.argv[xLabelField]

if '-ylabel' in sys.argv:
    yLabelField = sys.argv.index('-ylabel') + 1
    yLabel = sys.argv[yLabelField]

if '-bins' in sys.argv:
    binField = sys.argv.index('-bins') + 1
    bins = int(sys.argv[binField])


if '-field' in sys.argv:
    binField = sys.argv.index('-field') + 1
    binnedField = int(sys.argv[binField])

pointOffset = 0.
if '-xoffset' in sys.argv:
    pointOffset = float(sys.argv[sys.argv.index('-xoffset') + 1])

if '-doLog' in sys.argv:
    doLog = True
    try:
        logField = sys.argv.index('-doLog') + 1
        logBase = int(sys.argv[logField])
    except:
        logBase = 10
    print 'taking log%d of x datapoints' % logBase
    xLabel = 'log' + str(logBase) + '(' + xLabel + ')'

if '-binSize' in sys.argv:
    lengthField = sys.argv.index('-binSize') + 1
    binLength = float(sys.argv[lengthField])

if '-ymax' in sys.argv:
    yField = sys.argv.index('-ymax') + 1
    maxY = int(sys.argv[yField])

if '-title' in sys.argv:
    figTitle = sys.argv[sys.argv.index('-title') + 1]

if '-figsize' in sys.argv:
    sizes = sys.argv[sys.argv.index('-figsize') + 1].strip().split(',')
    figure(figsize=(float(sizes[0]),float(sizes[1])))

binLabels = []
doLabels = False
if '-binlabels' in sys.argv:
    binLabels = sys.argv[sys.argv.index('-binlabels') + 1].strip().split(',')
    doLabels = True

barsLegend = []
if '-legend' in sys.argv:
    barsLegend = sys.argv[sys.argv.index('-legend') + 1].strip().split(',')
    
ind2 = arange(bins)

bars = []
barsColors = []
index = 0
for aFile in fileList:
    distbin = bins * [0]
    
    dataList = []
    for line in aFile:
        fields = line.strip().split()
        try:
            point = float(fields[binnedField]) + pointOffset
            if doLog:
                if point < 1:
                    point = 1
                point = log(point, logBase)
            dataList.append(point)
        except:
            continue
    
    print '%d data points' % len(dataList)
    
    dataList.sort()
    print 'low = %f high = %f' % (dataList[0], dataList[-1])
    
    if binLength < 0:
        binLength = abs(dataList[-1] - dataList[0]) / bins
    for point in dataList:
        try:
            distbin[int(round(point/binLength))] += 1
        except:
            #print point, binLength, int(round(point/binLength))
            distbin[-1] += 1
    
    print binLength, int(round(point/binLength))

    bars.append(bar(ind2 + offset[index], distbin, width, color=colorList[index]))
    barsColors.append(bars[-1][0])
    
    print distbin
    halfCount = sum(distbin) / 2
    median = 0
    foundMedian = False
    while not foundMedian:
        if sum(distbin[:median]) < halfCount:
            median += 1
        else:
            foundMedian = True
    print median
    index += 1

xlim(-1 * width - 0.2, bins + 0.2)

if len(barsLegend) > 0:
    legend(barsColors, barsLegend)

ylabel(yLabel)
xlabel(xLabel)

if doLabels:
    setp(gca(), 'xticklabels', binLabels)

if maxY > 0:
    ylim(0, maxY)

if len(figTitle) > 0:
    title(figTitle)
gca().get_xaxis().tick_bottom()
gca().get_yaxis().tick_left()

savefig(pngfilename)
