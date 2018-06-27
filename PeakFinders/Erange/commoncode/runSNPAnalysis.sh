#!/bin/bash
#
# runSNPAnalysis.sh
#
# Usages: $ERANGEPATH/runSNPAnalysis.sh mouse rdsfile label rmaskdbfile dbsnpfile uniqStartMin totalRatio rpkmfile cachepages
# Example: /getSNPs.sh mouse /woldlab/trog/sdc/alim/24T4spike_10212/24T4spike.rds 24Tspike /woldlab/trog/data1/wlee/db/rmask.db /woldlab/trog/data1/wlee/db/dbSNP128.db 5 0.75 ~/proj/c2c12rna24R/c2c12rna.24R.final.rpkm 5000000 

# set ERANGEPATH to the absolute or relative path to ERANGE, if it's not in the environment

if [ -z "$ERANGEPATH" ]
then
    ERANGEPATH='../commoncode'
fi

echo 'runSNPAnalysis.sh: version 3.1'

cachepages=""
if [ $# -eq 9 ]; then
    cachepages="-cache "$9
fi

nosplices=""
if [ $# -eq 10 ]; then
    nosplices=" -nosplices "
fi

if [ $# -lt 8 ]; then
    echo 'runSNPAnalysis.sh genome rdsfile label rmaskdbfile dbsnpfile uniqStartMin totalRatio rpkmfile [cachepages]'
    echo 'where for each position S:'
    echo '     uniqStartMin = # independent reads supporting base change at S'
    echo '     totalRatio = total # reads supporting base change at S / total # reads that pass through S'
else
# log the parameters
arguments=$1' '$2' '$3' '$4' '$5' '$6' '$7' '$8' '$cachepages$nosplices
echo 'running with settings: ' $arguments
python $ERANGEPATH/recordLog.py snp.log runSNPAnalysis.sh "with parameters: $arguments"

# get all SNPs by extracting it from the RDS
python $ERANGEPATH/getSNPs.py $2 $6 $7 $3.snps.txt -enforceChr $cachepages $nosplices

# get SNPs in non-repeat regions only
python $ERANGEPATH/chkSNPrmask.py $4 $3.snps.txt $3.nr_snps.txt $cachepages

# Check to see if SNPs are found in dbSNP
# if dbSNP128.db is not built yet, build it by running buildsnpdb.py - build snp database using the dbSNP database file downloaded from UCSC
# usage: python2.5 buildsnpdb.py snpdbdir snpdbname
# the database flat file must be in the snpdbdir directory
# To build dbSNP database file, run the following command 
# python2.5 buildsnpdb.py snp128.txt dbSNP128

# get dbSNP info for SNPs that are found in the dbSNP database
python $ERANGEPATH/chksnp.py $5 $3.nr_snps.txt $3.nr_dbsnp.txt $cachepages

# get gene info for the snps found in dbSNP
python $ERANGEPATH/getSNPGeneInfo.py $1 $3.nr_dbsnp.txt $8 $3.nr_dbsnp_geneinfo.txt $cachepages

# get gene info for snps that are not found in dbSNP
python $ERANGEPATH/getNovelSNPs.py $1 $3.nr_dbsnp_geneinfo.txt $3.nr.final.txt 

# make bed file for displaying the snps on UCSC genome browser
python $ERANGEPATH/makeSNPtrack.py $3.nr_snps.txt $3 $3.nr_snps.bed
fi