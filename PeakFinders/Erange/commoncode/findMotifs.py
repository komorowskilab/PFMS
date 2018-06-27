#
#  findMotifs.py
#  ENRAGE
#
try:
    import psyco
    psyco.full()
except:
    pass

import sys, os
print '%s: version 3.4' % sys.argv[0]

numMotifs = '10'
maxWidth = 28
threshold = 75.

if len(sys.argv) < 3:
    print 'usage: python %s explabel regions.fsa [-meme] [-cisGreedy] [-logo] [-threshold percent] [-prefix motifprefix] [-numMotifs num] [-maxWidth bp] [-maskLower]' % sys.argv[0]
    print '\n\twhere at least one of the motif finders (meme or cisGreedy) must be specified\n'
    sys.exit(1)

expbase = sys.argv[1]
fsafile = sys.argv[2]

motifPrefix = expbase

doMeme = False
doCisGreedy = False
doCisSampler = False
noFinder = True
saveLogo = False
maskLower = False

if '-prefix' in sys.argv:
    motifPrefix = sys.argv[sys.argv.index('-prefix') + 1]

if '-logo' in sys.argv:
    saveLogo = True
    
if '-numMotifs' in sys.argv:
    numMotifs = sys.argv[sys.argv.index('-numMotifs') + 1]

if '-maxWidth' in sys.argv:
    maxWidth = int(sys.argv[sys.argv.index('-maxWidth') + 1])

if '-threshold' in sys.argv:
    threshold = float(sys.argv[sys.argv.index('-threshold') + 1])

if '-meme' in sys.argv:
    doMeme = True
    noFinder = False
    
if '-cisGreedy' in sys.argv:
    doCisGreedy = True
    noFinder = False

if '-cisSampler' in sys.argv:
    print 'cisSampler is not supported yet! avoid using it for now'
    doCisSampler = True
    noFinder = False

if '-maskLower' in sys.argv:
    maskLower = True
    
if noFinder:
    print 'error: must specify at least one motif finder - exiting'
    sys.exit(1)

from cistematic.experiments.fasta import Fasta

from cistematic.programs.meme import Meme
from cistematic.programs.cisGreedy import CisGreedy
#from cistematic.programs.alignace import AlignACE
from cistematic.programs.locator import Locator

exp = Fasta(expbase, expbase + '.db')
#exp.debugSQL = True

exp.initialize()
if maskLower:
    exp.setMaskLowerCase(True)

if doMeme:
    prog4 = Meme()
    prog4.setMaxWidth(maxWidth)
    prog4.setNumMotifs(numMotifs)
    prog4.setModel('zoops')
    exp.appendProgram(prog4)

if doCisGreedy:
    prog5 = CisGreedy()
    prog5.setGenExpOptions([])
    prog5.setMaxWidth(maxWidth)
    prog5.setNumMotifs(numMotifs)
    exp.appendProgram(prog5)

# cisSampler is not supported yet!
if doCisSampler:
    from cistematic.programs.cisSampler import CisSampler
    prog6 = CisSampler()
    prog6.setGenExpOptions([])
    prog6.setMaxWidth(maxWidth)
    prog6.setNumMotifs(numMotifs)
    exp.appendProgram(prog6)

exp.run(fsafile)

exp.createAnalysis()
exp.loadAnalysis()
exp.mapMotifs(threshold, verbose=False)
exp.exportMotifs(prefix = motifPrefix)
if saveLogo:
    exp.exportLogos(prefix = motifPrefix)

exp.draw(expbase + '.png', maxOccurences=4000)
print 'deleting database...'
del exp
os.remove(expbase + '.db')

