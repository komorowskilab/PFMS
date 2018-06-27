#
#  scatterfields.py
#  ENRAGE
#
import matplotlib
matplotlib.use('Agg')

from pylab import *
import math, cmath
import sys

doLogF1 = False
doLogF2 = False
doArcsinh = False
alphaVal = 0.5

print "%s: version 3.1" % sys.argv[0]
if len(sys.argv) < 7:
    print "usage: python %s infilename xaxisLabel xField yaxisLabel yField outImageName [-xmin xMin] [-ymin yMin] [-xmax xMax] [-ymax yMax] [-doLogF1] [-doLogF2] [-arcsinh] [-order polyOrder] [-base logBase] [-markGenes geneFile] [-markfold times] [-noregression] [-large] [-markdiag] [-title text] [-verbose]" % sys.argv[0]
    print "\n\tDo a scatter plot of 2 fields from an input file."
    print "\tfields are counted from 0."
    print "\tuse [-order polyOrder] to specify polynomial fits > 1"
    print "\tSupports very rudimentary compound fields for X value," 
    print "\tusing python's lambda functions (omit the keyword lambda)\n"
    sys.exit(1)

base = 10
fitOrder = 1
infile = open(sys.argv[1])
xaxis = sys.argv[2]
compoundField = False
doRegression = True
xField = sys.argv[3]
try:
    xField = int(xField)
except:
    try:
        compoundOp = 'lambda ' + xField
        operator = eval(compoundOp)
        compoundField = True
        print "compound field %s" % xField
    except:
        pass
    if not compoundField:
        print "expression %s not supported" % xField
        sys.exit(1)
yaxis = sys.argv[4]
yField = int(sys.argv[5])
outfilename = sys.argv[6]

forcexmax = -1
forceymax = -1

forcexmin = 0.
forceymin = 0.
figtitle = ''

if '-xmin' in sys.argv:
    xminField = sys.argv.index('-xmin') + 1
    forcexmin = float(sys.argv[xminField])

if '-ymin' in sys.argv:
    yminField = sys.argv.index('-ymin') + 1
    forceymin = float(sys.argv[yminField])

if '-xmax' in sys.argv:
    xmaxField = sys.argv.index('-xmax') + 1
    forcexmax = float(sys.argv[xmaxField])

if '-ymax' in sys.argv:
    ymaxField = sys.argv.index('-ymax') + 1
    forceymax = float(sys.argv[ymaxField])

if '-doLogF1' in sys.argv:
    doLogF1 = True

if '-doLogF2' in sys.argv:
    doLogF2 = True

if '-arcsinh' in sys.argv:
    doArcsinh = True
    
if '-order' in sys.argv:
    orderField = sys.argv.index('-order') + 1
    fitOrder = int(sys.argv[orderField])

if '-base' in sys.argv:
    baseField = sys.argv.index('-base') + 1
    base = int(sys.argv[baseField])

if '-title' in sys.argv:
    figtitle = sys.argv[sys.argv.index('-title') + 1]

verbose = False
if '-verbose' in sys.argv:
    verbose = True
    
markedGenes = []
marking = False
if '-markGenes' in  sys.argv:
    markField = sys.argv.index('-markGenes') + 1
    markFile = open(sys.argv[markField])
    for line in markFile:
        try:
            markedGenes.append(line.strip().split()[0].upper())
        except:
            markedGenes.append(line.strip().upper())        
    markFile.close()
    marking = True

markFold = False
folds = 1.
if '-markfold' in sys.argv:
    markFold = True
    foldChange = float(sys.argv[sys.argv.index('-markfold') + 1])
    
if '-noregression' in sys.argv:
    doRegression = False

plotLarge = False
if '-large' in sys.argv:
    plotLarge = True

markDiag = False
if '-markdiag' in sys.argv:
    markDiag = True
    print 'marking diagonal'

newscores = []
oldscores = []

markednewscores = []
markedoldscores = []

markedfoldnewscores = []
markedfoldoldscores = []

ymax = 0.
xmax = 0.
for line in infile:
    fields = line.strip().split()
    gene = fields[0]
    try:
        if compoundField:
            score = operator(fields)
        else:
            score = float(fields[xField])
        newscore = float(fields[yField])
    except:
        continue
    foldMarkThisScore = False
    if markFold:
        tempscore = score
        if tempscore == 0:
            tempscore = 0.03
        tempratio = newscore / tempscore
        if tempratio == 0:
            tempratio2 = tempscore / 0.03
        else:
            tempratio2 = 1. / tempratio
        if tempratio > foldChange or tempratio2 > foldChange:
            foldMarkThisScore = True
    if doArcsinh:
        score = abs(cmath.asinh(score))
    elif doLogF1:
        #if score < 1.0:
        #    score = 1.0
        #    #continue
        try:
            score = math.log(score, base)
        except:
            score = forcexmin
        if score > xmax:
            xmax = score
    if doArcsinh:
        newscore = abs(cmath.asinh(newscore))
    elif doLogF2:
        #newscore += 1
        #if newscore < 1.0:
        #    newscore = 1.0
        #    #continue
        try:
            newscore = math.log(newscore, base)
        except:
            newscore = forceymin
        if newscore > ymax:
            ymax = newscore
        #if xmax > 0. and newscore > xmax:
        #    xmax = newscore
    oldscores.append(score)
    newscores.append(newscore)
    if foldMarkThisScore:
        markedfoldoldscores.append(score)
        markedfoldnewscores.append(newscore)
        if marking and gene.upper() not in markedGenes:
            print gene, score, newscore, 'unmarked'
        if gene.upper() in markedGenes:
            print gene, score, newscore, 'overfold'
        if verbose:
            print len(markedfoldoldscores), line.strip()
    if gene.upper() in markedGenes:
        if not foldMarkThisScore:
            print gene, score, newscore
        markedoldscores.append(score)
        markednewscores.append(newscore)
    #print score, newscore
print score, newscore
print fields

if plotLarge and markFold:
    plot(oldscores, newscores, '^', markersize=10., color='0.75', alpha=alphaVal)
elif plotLarge:
    plot(oldscores, newscores, 'b^', markersize=10., alpha=alphaVal)
elif markFold:
    plot(oldscores, newscores, ',', color='0.75', alpha=alphaVal)
else:
    plot(oldscores, newscores, 'b,', alpha=alphaVal)

if len(markedfoldoldscores) > 0:
    if plotLarge:
        plot(markedfoldoldscores, markedfoldnewscores, 'b^', markersize=10., alpha=alphaVal)
    else:
        plot(markedfoldoldscores, markedfoldnewscores, 'b,', alpha=alphaVal)

if len(markedoldscores) > 0:
    if plotLarge:
        plot(markedoldscores, markednewscores, 'r^', color='red', markersize=10., alpha=alphaVal)
    else:
        plot(markedoldscores, markednewscores, '.', color='red', markersize=4., alpha=alphaVal)

fitvalues = polyfit(oldscores, newscores, fitOrder)
print fitvalues
#fitvalues = polyfit(oldscores, newscores, 1)
print len(oldscores)

meanObserved = float(sum(newscores)) / len(newscores)
if len(fitvalues) == 2:
    predicted = [(fitvalues[0] * x + fitvalues[1]) for x in oldscores]
else:
    predicted = [(fitvalues[0] * x**2 + fitvalues[1] * x + fitvalues[2]) for x in oldscores]

meanPredicted = float(sum(predicted)) / len(predicted)
SSt = 0.
SSe = 0.

for index in range(len(newscores)):
    SSt += (newscores[index] - meanObserved) ** 2
    SSe += (newscores[index] - predicted[index]) ** 2

rSquared = 1. - SSe / SSt
print 'R**2 = %f' % rSquared
    
oldscores.sort()
if len(fitvalues) == 2:
    predicted = [(fitvalues[0] * x + fitvalues[1]) for x in oldscores]
else:
    predicted = [(fitvalues[0] * x**2 + fitvalues[1] * x + fitvalues[2]) for x in oldscores]
if doRegression:
    plot(oldscores, predicted, '-k', linewidth=2)

if figtitle == '':
    figtitle = '%s vs %s (R^2: %.2f)' % (yaxis, xaxis, rSquared)
title(figtitle)

if markDiag:
    min = forcexmin
    if forceymin < min:
        min = forceymin
    max = xmax
    if ymax > max:
        max = ymax
    if forcexmax > max:
        max = forcexmax
    if forceymax > max:
        max = forceymax
    plot([min,max], [min,max], '-g', linewidth=2)

print forcexmin, forceymin

if doLogF2:
    ylabel('log' + str(base) + '(' + yaxis + ')')
else:
    ylabel(yaxis)

if doLogF1:
    xlabel('log' + str(base) + '(' + xaxis + ')')
else:
    xlabel(xaxis)

if xmax > 0:
    xlim(forcexmin - 0.05,xmax)
if ymax > 0:
    ylim(forceymin - 0.05,ymax)
if forcexmax > 0 and forceymax > 0:
    xlim(forcexmin - 0.05, forcexmax)
    ylim(forceymin - 0.05, forceymax)

gca().get_xaxis().tick_bottom()
gca().get_yaxis().tick_left()

savefig(outfilename, dpi=100)
