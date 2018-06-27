#
#  stallCategory.py
#  ENRAGE
#

try:
        import psyco
        psyco.full()
except:
        pass

import sys

print '%s: version 1.1' % sys.argv[0]

if len(sys.argv) < 4:
    print 'usage: python %s stalledPercentFile1 stalledPercentFile2 transcriptFile [-out oufile] [-statout statoutfile] [-expression level]' % sys.argv[0]
    sys.exit(1)

infile1 = open(sys.argv[1])
infile2 = open(sys.argv[2])
transcriptFile = open(sys.argv[3])

writeOut = False
if '-out' in sys.argv:
    outfile = open(sys.argv[sys.argv.index('-out') + 1],'w')
    outfile.write('gene\texpression\tratio1\tpromAmount1\ttotal1\trestRPKM1\tratio2\tpromAmount2\ttotal2\trestRPKM2\n')
    writeOut = True

statWriteOut = False
if '-statout' in sys.argv:
    statoutfile = open(sys.argv[sys.argv.index('-statout') + 1],'w')
    statoutfile.write('ExpressionR1R2Stalled1Stalled2\tCount\n')
    statWriteOut = True

expressionLevel = 0.9
if '-expression' in sys.argv:
    expressionLevel = float(sys.argv[sys.argv.index('-expression') + 1])

dictOne = {}
dictTwo = {}
expressionDict = {}

for line in infile1:
    if 'short' in line:
        continue
    fields = line.strip().split()
    promAmount = float(fields[4]) + float(fields[5])
    genelen = float(fields[3])/100
    total = float(fields[2])
    if total < 0.1:
        total = 0.1
    restRPKM = (total * (1. - promAmount/100.))/ (genelen - 0.6)
    
    ratio = float(fields[-1])
    dictOne[fields[1]] = (ratio, promAmount, total, restRPKM)

for line in infile2:
    if 'short' in line:
        continue
    fields = line.strip().split()
    promAmount = float(fields[4]) + float(fields[5])
    genelen = float(fields[3])/100
    if promAmount == 0.:
        promAmount = 0.1
    total = float(fields[2])
    if total < 0.1:
        total = 0.1
    restRPKM = (total * (1. - promAmount/100.))/ (genelen - 0.6)
    ratio = float(fields[-1])
    dictTwo[fields[1]] = (ratio, promAmount, total, restRPKM)

for line in transcriptFile:
    (gene, transc, transcpercell) = line.strip().split()
    expressionDict[gene] = float(transcpercell)

categoryList = []
categoryDict = {}
for atype in ['HH', 'HL', 'LH', 'LL']:
    for expression in ['E','N']:
        for cat1 in ['Y', 'N']:
            for cat2 in ['Y', 'N']:
                category = expression + cat1 + cat2 + atype
                categoryList.append(category)
                categoryDict[category] = []
for gene in dictOne:
    if gene not in expressionDict:
        if writeOut:
            print '%s is not in expressionDict - skipping' % gene
        continue
    expression = expressionDict[gene]
    (ratio1, promAmount1, total1, restRPKM1) = dictOne[gene]
    (ratio2, promAmount2, total2, restRPKM2) = dictTwo[gene]
    category = 'test'
    
    #
    if expression > expressionLevel:
        category = 'E'
    else:
        category = 'N'
    
    if total1 > 5.0:
        category += 'Y'
    else:
        category += 'N'
    
    if total2 > 5.0:
        category += 'Y'
    else:
        category += 'N'
    
    if ratio1 > 15:
        category += 'H'
    else:
        category += 'L'
    
    if ratio2 > 15:
        category += 'H'
    else:
        category += 'L'
    
    categoryDict[category].append(gene)
    if writeOut:
        outfile.write('%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n' % (gene, expression, ratio1, promAmount1, total1, restRPKM1, ratio2, promAmount2, total2, restRPKM2, category)
)
    else:
        print '%s %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %s' % (gene, expression, ratio1, promAmount1, total1, restRPKM1, ratio2, promAmount2, total2, restRPKM2, category)

if writeOut:
    outfile.close()

for category in categoryList:
    if statWriteOut:
        statoutfile.write('%s\t%d\n' % (category, len(categoryDict[category])))
    else:
        print '%s %d' % (category, len(categoryDict[category]))


