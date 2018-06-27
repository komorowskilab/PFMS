import sys

print "%s: version 1.0" % sys.argv[0]
if len(sys.argv) < 3:
    print "usage: python %s infile.gff outfile.cis\n" % sys.argv[0]
    print "\tTHIS SCRIPT WILL MOST LIKELY NEED TO BE EDITED FOR YOUR GFF FILE\n"
    sys.exit(1)

index = 1
# Cistematic just want's a use set of exons labeled "CDS", "5UTR", and "3UTR"
# just put the corresponding type in your GFF file as the key in the key:value pairs 
# in the ftypeDict below
ftypeDict = {"CDS":"CDS", "mRNA":"mRNA", "five_prime_utr":"5UTR", "three_prime_utr":"3UTR"}
chrom = ''
idfields = ''
gene = ''
sense = ''
start = 0
stop = 0
ftype = ''

infile = open(sys.argv[1])
outfile = open(sys.argv[2],'w')
for line in infile:
    if line[0]=='#':
        continue
    fields = line.strip().split()
    try:
        if fields[2] in ftypeDict:
            # this part of the code will need to be customized, most likely
            # how does the annotation define the gene, geneid, and chromosome
            # for example, for Anopheles Gambiae we have
            #chrX    VectorBase      mRNA    582     16387   .       -       .       ID=vectorbase|AGAP000002-RA; stable_id=AGAP000002-RA.1; Parent=vectorbase|AGAP000002;
            if fields[2] == 'mRNA':
                #print fields
                chrom = fields[0][3:]
                source = fields[1]
                idfields = fields[9].split(';')
                geneid = idfields[0].split('=')[1]
                sense = fields[6]
                #print chrom, idfields, gene, sense
            else:
                start = int(fields[3])
                stop = int(fields[4])
                ftype = ftypeDict[fields[2]]
                outline = '%s\t%s%d\t%s\t%d\t%d\t%s\t%s\n' % (geneid, source, index, chrom, start, stop, sense, ftype)
                #outline = '%s\tAUGUSTUS%d\tchr%s\t%d\t%d\t%s\t%s\n' % (gene, index, chrom, start, stop, sense, ftype)
                outfile.write(outline)
    except:
        sys.exit()
    index += 1
infile.close()
outfile.close()
