#
#  plotnomogram.py
#  ENRAGE
#

import sys
print '%s: version 1.1' % sys.argv[0]
import matplotlib
matplotlib.use('Agg')

from pylab import *
import matplotlib.axes 

try:
	import psyco
	psyco.full()
except:
	pass

if len(sys.argv) < 5:
    print 'usage: python %s maxdev xreads infile outpng' % sys.argv[0]
    sys.exit(1)

maxdev = float(sys.argv[1])
xreads = float(sys.argv[2])
infile = open(sys.argv[3])
outfilename = sys.argv[4]

line = infile.readline().strip()

percentages = line.split()
del percentages[0]

listWidth = len(percentages)

geneValues = {}

for line in infile:
    fields = line.strip().split()
    geneValues[fields[0]] = []
    for pos in range(listWidth):
        geneValues[fields[0]].append(float(fields[1 + pos]))
        
# categories here are: 3000+, 2999-300, 299-30, 29-3
genes3000p = []
genes300p = []
genes30p = []
genes3p = []

for gene in geneValues:
    finalLevel = geneValues[gene][0]
    if finalLevel >= 3000:
        genes3000p.append(gene)
    elif finalLevel >= 300:
        genes300p.append(gene)
    elif finalLevel >= 30:
        genes30p.append(gene)
    elif finalLevel >= 3:
        genes3p.append(gene)

organizedList = [genes3000p, genes300p, genes30p, genes3p]
listNames = ['3000+ RPKM     ', '300-2999 RPKM', '30-299 RPKM    ', '3-29 RPKM        ']
listColors = ['k','c', 'm', 'r']
geneCounts = {}
oldscores = [0.]
newscores = {}
for name in listNames:
    newscores[name] = [0.]

index = 0
for percent in percentages[1:]:
    oldscores.append(xreads * float(percent) / 100.)
    index += 1
    listindex = 0
    for geneList in organizedList:
        geneCount = len(geneList)
        numOver = 0.
        for gene in geneList:
            finalVal = geneValues[gene][0]
            currentVal = geneValues[gene][index]
            if abs((currentVal - finalVal) / finalVal) > maxdev:
                numOver += 1.
        fraction = 1. - numOver / geneCount
        print '%s %s %d %.2f' % (percent, listNames[listindex], geneCount, fraction)
        newscores[listNames[listindex]].append(fraction)
        geneCounts[listNames[listindex]] = geneCount        
        listindex += 1

matplotlib.axes._process_plot_var_args.defaultColors = ['k','y','m','c','b','g','r']

oldscores.append(xreads)
index = 0
plots = []
plotsColors = []
plotsLegend = []
for name in listNames:
    newscores[name].append(1.0)
    plots.append(plot(oldscores, newscores[name], listColors[index], linewidth=2))
    plot(oldscores[1:-1], newscores[name][1:-1], listColors[index] + '^')
    plotsColors.append(plots[-1][0])
    plotsLegend.append('%s n = %d' % (name, geneCounts[name]))
    index += 1
legend(plotsColors, plotsLegend, loc=0)
xticks(oldscores)
locs, labels = xticks()
setp(labels, rotation='vertical')
ylim(0, 1.03)
xlim(-0.1, xreads + .1)
savefig(outfilename)
