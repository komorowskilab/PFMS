#
#  profilebins.py
#  ENRAGE
#

try:
    import psyco
    psyco.full()
except:
    pass

import sys
print '%s: version 2.2' % sys.argv[0]

from commoncode import *

if len(sys.argv) < 4:
    print 'usage: python %s label infile1 [-upstream infile2] [-downstream infile3] [-uplength kb] [-downlength kb] [-gene geneName] [-genes genefile] [-append] outfile' % sys.argv[0]
    sys.exit(1)

label = sys.argv[1]
infilename1 = sys.argv[2]
fileList = [infilename1]

outfilename = sys.argv[-1]

geneList = []

restrictGenes = False

if '-gene' in sys.argv:
    geneList.append(sys.argv[sys.argv.index('-gene') + 1])
    restrictGenes = True
if '-genes' in sys.argv:
    geneindex = sys.argv.index('-genes') + 1
    genefile = open(sys.argv[geneindex])
    for line in genefile:
        fields = line.strip().split()
        if len(fields) > 1:
            geneList.append(fields[0])
        else:
            geneList.append(line.strip())
    restrictGenes = True

doUp = False
if '-upstream' in sys.argv:
    upfilefield = sys.argv.index('-upstream') + 1
    upfilename = sys.argv[upfilefield]
    doUp = True
    fileList = [upfilename, infilename1]

doDown = False
if '-downstream' in sys.argv:
    downfilefield = sys.argv.index('-downstream') + 1
    downfilename = sys.argv[downfilefield]
    doDown = True
    fileList.append(downfilename)

doAppend = False
if '-append' in sys.argv:
    doAppend = True

partLength = [10.]
partOffset = [0.]

if '-uplength' in sys.argv:
    uplengthField = sys.argv.index('-uplength') + 1
    uplength = float(sys.argv[uplengthField])
    partLength = [uplength, 10.]
    partOffset = [-1. * uplength, 0.]

if '-downlength' in sys.argv:
    downlengthField = sys.argv.index('-downlength') + 1
    downlength = float(sys.argv[downlengthField])
    partLength.append(downlength)
    partOffset.append(10.)

totalWeight = 0.
totalBins = []
for afile in fileList:   
    infile = open(afile)
    
    line = infile.readline()
    fields = line.strip().split()
    numBins = len(fields) - 4

    geneName = fields[1]
    weight = float(fields[2])
    if restrictGenes and geneName in geneList:
        totalWeight += weight
    totalBins.append([])
    for myBin in fields[4:]:
        if not restrictGenes or (restrictGenes and geneName in geneList):
            totalBins[-1].append(weight * float(myBin))
        else:
            totalBins[-1].append(0.)
    geneDict = {}
    
    for line in infile:
        fields = line.strip().split()
        geneName = fields[1]
        if restrictGenes and geneName not in geneList:
            continue
        weight = float(fields[2])
        index = 0
        for myBin in fields[4:]:
            totalBins[-1][index] += weight * float(myBin)
            index += 1
        totalWeight += weight

sumWeight = 0.
totalPercent = 0.
if doAppend:
    outfile = open(outfilename,'a')
else:
    outfile = open(outfilename,'w')
    outfile.write('x-axis')
    partIndex = 0
    for partBins in totalBins:
        partLen = partLength[partIndex]
        numBins = len(partBins)
        for binIndex in range(numBins):
            outfile.write('\t%.2f' % (partOffset[partIndex] + (binIndex * partLen/numBins)))
        partIndex += 1
    outfile.write('\tweight\n')

outfile.write(label)
for partBins in totalBins:
    for aBin in partBins:
        percent = aBin / totalWeight
        outfile.write('\t%.1f' % percent)
        sumWeight += aBin
        totalPercent += percent
outfile.write('\t%.1f\n' % totalWeight)
outfile.close()

print sumWeight
print totalPercent
