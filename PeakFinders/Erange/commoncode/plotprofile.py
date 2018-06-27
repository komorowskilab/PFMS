#
#  plotprofile.py
#  ENRAGE
#

import sys
try:
    import psyco
    psyco.full()
except:
	pass

print '%s: version 2.2' % (sys.argv[0])
if len(sys.argv) < 3:
    print 'usage: python %s infile outfile.png [-scale] [-max weightMax] [-ymin bottom] [-ymax top] [-subtractEvens]' % sys.argv[0]
    sys.exit(1)

doScale = False
weightMax = 1.

if '-scale' in sys.argv:
    doScale = True

if '-max' in sys.argv:
    weightMax = float(sys.argv[sys.argv.index('-max') + 1])

limitYscale = False
ymax = 0.
if '-ymax' in sys.argv:
    ymax = float(sys.argv[sys.argv.index('-ymax') + 1])
    limitYscale = True

ymin = 0.
if '-ymin' in sys.argv:
    ymin = float(sys.argv[sys.argv.index('-ymin') + 1])
    limitYscale = True

subtractEvens = False
if '-subtractEvens' in sys.argv:
    subtractEvens = True

import matplotlib
matplotlib.use('Agg')

from pylab import *
from math import *

infile = open(sys.argv[1])
pngfilename = sys.argv[2]

labelList = []
dataList = []
plotList = []
xmin = 10**20
xmax = -10**20

xcoordList = []
datapointList = []
weightList = []
line = infile.readline()
fields = line.strip().split()
for data in fields[1:-1]:
    datapoint = float(data)
    if datapoint < xmin:
        xmin = datapoint
    if datapoint > xmax:
        xmax = datapoint
    xcoordList.append(datapoint)

index = 1
for line in infile:
    fields = line.strip().split()
    datapointList = []
    for data in fields[1:-1]:
        datapointList.append(float(data))
    if subtractEvens and index % 2 == 0:
        for dataIndex in range(len(datapointList)):
            dataList[-1][dataIndex] -= datapointList[dataIndex]
    else:
        dataList.append(datapointList)
    weight = float(fields[-1])
    if subtractEvens and index % 2 == 0:
        pass
    else:
        labelList.append(fields[0])
        if weight > weightMax:
            weightMax = weight
        weightList.append(weight)
    index += 1
    
for index in range(len(dataList)):
    newList = []
    if doScale:
        scale = weightList[index] / weightMax
        print weightList[index], weightMax, scale
        for val in dataList[index]:
            newList.append(val * scale)
    else:
        newList = dataList[index]
    plotList.append(plot(xcoordList, newList, linewidth=3.0))

xticks(xcoordList, rotation='vertical')
xlim(xmin - 0.1, xmax + 0.1)
if limitYscale:
    ylim(ymin, ymax)

legend(plotList, labelList)
savefig(pngfilename)
