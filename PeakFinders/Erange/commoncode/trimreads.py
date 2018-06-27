#
#  trimquery.py
#  ENRAGE
#
#  Created by Ali Mortazavi on 8/12/08.
#

import sys
from cistematic.core import complement

print '%s: version 2.1' % sys.argv[0]
if len(sys.argv) < 4:
    print "usage: python %s length infile outfile [-fastq] [-fromback] [-paired] [-flip] [-filter maxN]" % sys.argv[0]
    print "\t where paired fragments are separated by a : when given the -paired flag" 
    sys.exit(1)

length = int(sys.argv[1])
infile = open(sys.argv[2])
outfile = open(sys.argv[3],'w')

fastq = False
if '-fastq' in sys.argv:
    fastq = True
    
paired = False
if '-paired' in sys.argv:
    paired = True
    pairedlength = 2 * length
index = 0
fromBack = False
if '-fromback' in sys.argv:
    fromBack = True
    length = -1 * length

filtering = False
maxN = 2
if '-filter' in sys.argv:
    filtering = True
    try:
        maxN = int(sys.argv[sys.argv.index('-filter') + 1])
    except:
        maxN = 2
    print "filtering out reads with more than %d Ns" % maxN

flipseq = False
if '-flip' in sys.argv:
    flipseq = True

print "trimming reads from %s to %d bp and saving them in %s" % (sys.argv[2], length, sys.argv[3])
    
filtered = 0
header = ''
for line in infile:
    line = line.strip()
    if len(line) == 0:
        continue
    firstChar = line[0]
    if (not fastq and firstChar == '>') or (fastq and firstChar in ['@','+']): 
        header = line + '\n'
    else:
        if filtering:
            if line.count('N') > maxN:
                filtered += 1
                continue
        seq1 = line[length:]
        seq2 = line[:length]
        if flipseq:
            try:
                tempseq1 = seq1
                seq1 = complement(tempseq1)
            except:
                seq1 = tempseq1
            try:
                tempseq2 = seq2
                seq2 = complement(tempseq2)
            except:
                seq2 = tempseq2
        if paired:
            if len(line) < pairedlength:
                continue
            outfile.write(header + seq1 + ':' + seq2 + '\n')
        else:
            if fromBack:
                outfile.write(header + seq1 + '\n')
            else:
                outfile.write(header + seq2 + '\n')
        index += 1
        if index % 1000000 == 0:
            print '.',
        sys.stdout.flush()
outfile.close()
print "returned %d reads" % index
if filtering:
    print "%d additional reads filtered" % filtered