try:
	import psyco
	psyco.full()
except:
	pass
from cistematic.genomes import Genome
from math import log
import os.path
import sys

print '%s: version 2.1' % sys.argv[0]

if len(sys.argv) < 6:
    print 'usage: python %s genome outimage gofileroot1 title1 cohortsize1 [gofileroot2 title2 cohortsize2 ...] [-fontsize pts] [-length in] [-width in]' % sys.argv[0]
    sys.exit(1)

hg = Genome(sys.argv[1])
allgodesc = hg.allGOterms()
godesc = []

import matplotlib
matplotlib.use('Agg')

from pylab import *
doGray = False

rootdir = './'

imagename =  sys.argv[2]

options = 0
fontSize = 5
length = 10
width = 7
if '-fontsize' in sys.argv:
    fontSize = int(sys.argv[sys.argv.index('-fontsize') + 1])
    options += 2

if '-length' in sys.argv:
    length = int(sys.argv[sys.argv.index('-length') + 1])
    options += 2

if '-width' in sys.argv:
    width = int(sys.argv[sys.argv.index('-width') + 1])
    options += 2

conditions = (len(sys.argv) - 3 - options) / 3
fileroots = []
titles = []
for index in range(conditions):
    fileroots.append(sys.argv[ 3 + index * 3])
    titles.append((sys.argv[4 + index * 3], '(' + sys.argv[5 + index * 3] + ')'))

htmlname = imagename[:-4] + '.html'

ceiling = 40.0
goterms = []
goscores = {}
numgenes = {}
possiblegenes = {}
flatArray = []

highestPval = 0.0
lowestPval = 1.0
for sigfile in fileroots:
    infile = open(rootdir + sigfile + '.gosig', 'r')
    for line in infile:
        if 'depleted' in line:
            continue
        fields = line.split('\t')
        if fields[0] not in goterms:
            goterms.append(fields[0])
            goscores[fields[0]] = []
            numgenes[fields[0]] = []
            possiblegenes[fields[0]] = 0
        if float(fields[3]) > highestPval:
            highestPval = float(fields[3])
        if float(fields[3]) < lowestPval:
            lowestPval = float(fields[3])

print highestPval
print lowestPval

boundaryScore = score = -1 * log(highestPval) /  (2.0 * ceiling) + 0.49
print boundaryScore

#cdict = {'red': ((0.0, 1.0, 1.0), (boundaryScore, 1.0, 0), (1.0, 1.0, 1.0)),
#	'green': ((0.0, 1.0, 1.0), (boundaryScore, 1.0, 0), (1.0, 1.0, 0.0)),
#	'blue': ((0.0, 1.0, 1.0), (boundaryScore, 1.0, 0), (1.0, 0.0, 0.0))}

cdict = {'red': ((0.0, 1.0, 1.0), (boundaryScore, 1.0, 0.1), (1.0, 1.0, 1.0)),
	'green': ((0.0, 1.0, 1.0), (boundaryScore, 1.0, 0.1), (1.0, 1.0, 1.0)),
	'blue': ((0.0, 1.0, 1.0), (boundaryScore, 1.0, 0.75), (1.0, 0.0, 0.0))}
mymap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,1024)

goindex = 0
for zfile in fileroots:
    infile = open(rootdir + zfile + '.gozscore', 'r')
    for line in infile:
        fields = line.split()
        goindex += 1
        if fields[0] not in goterms:
            continue
        score = -1 * log(float(fields[7])) /  (2.0 * ceiling)
        if score < -0.5:
            score = -0.5
        if score > 0.5:
            score = 0.5
        score += 0.5
        if doGray:
            score = 1 - score
        goscores[fields[0]].append(score)
        numgenes[fields[0]].append(fields[1])
        possiblegenes[fields[0]] = int(fields[4])
goindex /= len(fileroots)

gokeys = goscores.keys()
gosortarray = []
for term in gokeys:
    gosortarray.append(goscores[term] + [term])
gosortarray.sort()

htmlfile = open(htmlname, 'w')
htmlfile.write('<html><head><title>GO Analysis</title></head><body><table border="1">')
htmlfile.write('<tr><th>Description</th><th>possible</th>')
for entry in titles:
    htmlfile.write('<th>%s<br>%s</th>' % entry)
htmlfile.write('</tr>\n')
tableLines = []

for entry in gosortarray:
    term = entry[-1]
    outline = term + ':\t'
    for entry in goscores[term]:
        outline += str(round(entry, 4)) + '\t'
    print outline
    htmlLine = '<tr><th>%s</th><th>%d</th>' % (allgodesc[term], possiblegenes[term])
    index = 0
    for fileroot in fileroots:
        gofile = fileroot + '.' + term[3:]
        ngene = numgenes[term][index]
        if os.path.exists(gofile):
            htmlLine += '<td><a href="%s">%s</a></td>' % (gofile, ngene)
        else:
            htmlLine += '<td>%s</td>' % (ngene)
        index += 1
    tableLines.append(htmlLine + '</tr>\n')
    flatArray.append(goscores[term])
    godesc.append(allgodesc[term])

tableLines.reverse()
for line in tableLines:
    htmlfile.write(line)
htmlfile.write('<tr><th>Cohort Size:</th>')
htmlfile.write('</tr>\n')
htmlfile.write('</table></body></html>')

figure(figsize=(length, width))
myaxe = axes([0.3, 0.1, 0.55, 0.75])

Z = array(flatArray)
print Z.shape
if doGray:
    c = pcolor(Z, cmap=cm.gray, vmin=0.0, vmax=1.0)
else:
    c = pcolor(Z, cmap=mymap, vmin=0.0, vmax=1.0)
c.set_linewidth(0.1)
clim(0.0, 1.0)

ind = arange(len(fileroots))
width = 0.5

coordy = 0.1
#deltaX = 6.0 * 0.85 / float(len(fileroots))
#deltaY = 7.0  * 4.12 / len(gosortarray)
deltaX = 1.0	
deltaY = 1.0

pcolorAxes = c.get_axes()
for entry in gosortarray:
    term = entry[-1]
    coordx = 0.4
    for genenum in numgenes[term]:
        if len(genenum) == 1:
            genenum = '    ' + genenum
        elif len(genenum) == 2:
            genenum = '  ' + genenum
        pcolorAxes.text(coordx, coordy, genenum, fontsize=fontSize)
        coordx += deltaX
    coordy += deltaY

coordx = 0
for (line1,line2) in titles:
    pcolorAxes.text(coordx + 0.1, coordy + 3 * deltaY + 0.5, line1, fontsize=int(fontSize*1.5))
    pcolorAxes.text(coordx + 0.1, coordy + deltaY, line2, fontsize=int(fontSize*1.5))
    coordx += deltaX 

setp(gca(), 'xticks', [])
setp(gca(), 'xticklabels', [])
setp(gca(), 'yticks', arange(len(godesc)))
setp(gca(), 'yticklabels', godesc)
locs, labels = yticks()
setp(labels, fontsize=fontSize)
setp(labels, verticalalignment='bottom')
setp(gca(), 'ylim', [0, len(godesc)])

figtext(0.3,0.02, str(goindex - len(gokeys)) + ' additional GO Terms below threshold of significance', fontsize=fontSize*2)

#d = colorbar(orientation="vertical", edgecolor='w', drawedges=False)
d = colorbar(orientation="vertical", drawedges=False)
for t in d.ax.get_yticklabels():
    t.set_fontsize(0)
locs, labels = yticks()
setp(labels, fontsize=5)
pcolorAxes.text(conditions + 1,len(godesc), str(lowestPval), fontsize=fontSize*2)
pcolorAxes.text(conditions + 1,boundaryScore * len(godesc), str(highestPval), fontsize=fontSize*2)

savefig(imagename, dpi=250)
show()

