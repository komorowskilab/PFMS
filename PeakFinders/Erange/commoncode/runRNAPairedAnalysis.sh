#!/bin/bash
#
# runRNAPairedAnalysis.sh
# ENRAGE
#
# example: . ../commoncode/runRNAPairedAnalysis.sh mouse c2c12rna ../mm9repeats/rmask.db
#       
#          assuming that we have rds database with the prefix c2c12rna.24R and that an RNAFAR analysis has already been run. 

# set ERANGEPATH to the absolute or relative path to ERANGE, if it's not in the environment

if [ -z "$ERANGEPATH" ]
then
    ERANGEPATH='../commoncode'
fi

echo 'runRNAPairedAnalysis.sh: version 3.7'

models=""
if [ $# -eq 5 ]; then
    models=" -models "$5
fi

replacemodels=""
if [ $# -eq 6 ]; then
    replacemodels=" -models $5 -replacemodels "
fi

if [ -z "$1" ]
then
    echo
    echo 'usage:runRNAPairedAnalysis.sh genome rdsprefix repeatmaskdb [modelfile] [-replacemodels]'
    echo
    echo 'where rdsprefix is the name of the rds file without the .rds extension'
    echo 'use "none" for the repeatmaskdb if you do not have one'
    echo
else

# log the parameters
arguments=$1' '$2' '$3' '$models' '$5
echo 'running with settings: ' $arguments
python $ERANGEPATH/recordLog.py rna.log runRNAPairedAnalysis.sh "with parameters: $arguments"

# count the unique reads falling on the gene models ; the nomatch files are 
# mappable reads that fell outside of the Cistematic gene models and not the 
# unmappable of Eland (i.e, the "NM" reads)
echo "python $ERANGEPATH/geneMrnaCounts.py $1 $2.rds $2.uniqs.count -markGID -cache 1 $models $replacemodels"
python $ERANGEPATH/geneMrnaCounts.py $1 $2.rds $2.uniqs.count -markGID -cache 1 $models $replacemodels

# calculate a first-pass RPKM to re-weigh the unique reads,
# using 'none' for the splice count
python $ERANGEPATH/normalizeExpandedExonic.py $1 $2.rds $2.uniqs.count none $2.firstpass.rpkm -cache $models $replacemodels

# recount the unique reads with weights calculated during the first pass
python $ERANGEPATH/geneMrnaCountsWeighted.py $1 $2.rds $2.firstpass.rpkm $2.uniqs.recount -uniq -cache 1 $models $replacemodels

# count splice reads
python $ERANGEPATH/geneMrnaCounts.py $1 $2.rds $2.splices.count -splices -noUniqs -markGID -cache 1 $models $replacemodels

# find new regions outside of gene models with reads piled up 
python $ERANGEPATH/findall.py RNAFAR $2.rds $2.newregions.txt -RNA -minimum 1 -nomulti -flag NM -log rna.log -cache 1

# filter out new regions that overlap repeats more than a certain fraction
python $ERANGEPATH/checkrmask.py $3 $2.newregions.txt $2.newregions.repstatus $2.newregions.checked -startField 1 -log rna.log -cache 1

# calculate the read densities
python $ERANGEPATH/regionCounts.py $2.newregions.checked $2.rds $2.newregions.good -markRDS -cache -log rna.log

# map all candidate regions that have paired ends overlapping with known genes
python $ERANGEPATH/rnafarPairs.py $1 $2.newregions.good $2.rds $2.candidates.txt -cache $models $replacemodels

# calculate expanded exonic read density
python $ERANGEPATH/normalizeExpandedExonic.py $1 $2.rds $2.uniqs.recount $2.splices.count $2.expanded.rpkm $2.candidates.txt $2.accepted.rpkm -cache $models $replacemodels

# weigh multi-reads
python $ERANGEPATH/geneMrnaCountsWeighted.py $1 $2.rds $2.expanded.rpkm $2.multi.count -accept $2.accepted.rpkm -multi -cache 1 $models $replacemodels

# calculate final exonic read density
python $ERANGEPATH/normalizeFinalExonic.py $2.rds $2.expanded.rpkm $2.multi.count $2.final.rpkm -multifraction -withGID -cache 

fi
