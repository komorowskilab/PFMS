#!/bin/bash
#
# runStandardAnalysis.sh
# ENRAGE
#
# example: . $ERANGEPATH/runStandardAnalysis.sh mouse c2c12rna ../mm9repeats/rmask.db 20000
#       
#          assuming that we have rds database with the prefix c2c12rna.24R and that an RNAFAR analysis has already been run. 

# set ERANGEPATH to the absolute or relative path to ERANGE, if it's not in the environment

if [ -z "$ERANGEPATH" ]
then
    ERANGEPATH='../commoncode'
fi

echo 'runStandardAnalysis.sh: version 4.2'

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
    echo 'usage:runStandardAnalysis.sh genome rdsprefix repeatmaskdb bpradius [modelfile] [-replacemodels]'
    echo
    echo 'where rdsprefix is the name of the rds file without the .rds extension'
    echo 'use "none" for the repeatmaskdb if you do not have one'
    echo
else

# log the parameters
arguments=$1' '$2' '$3' '$4' '$models' '$6
echo 'running with settings: ' $arguments
python $ERANGEPATH/recordLog.py rna.log runStandardAnalysis.sh "with parameters: $arguments"

# count the unique reads falling on the gene models ; the nomatch files are 
# mappable reads that fell outside of the Cistematic gene models and not the 
# unmappable of Eland (i.e, the "NM" reads)
echo "python $ERANGEPATH/geneMrnaCounts.py $1 $2.rds $2.uniqs.count -markGID -cache 1 $models $replacemodels"
python $ERANGEPATH/geneMrnaCounts.py $1 $2.rds $2.uniqs.count -markGID -cache 1 $models $replacemodels

# calculate a first-pass RPKM to re-weigh the unique reads,
# using 'none' for the splice count
echo "python $ERANGEPATH/normalizeExpandedExonic.py $1 $2.rds $2.uniqs.count none $2.firstpass.rpkm -cache  $models $replacemodels"
python $ERANGEPATH/normalizeExpandedExonic.py $1 $2.rds $2.uniqs.count none $2.firstpass.rpkm -cache  $models $replacemodels

# recount the unique reads with weights calculated during the first pass
echo "python $ERANGEPATH/geneMrnaCountsWeighted.py $1 $2.rds $2.firstpass.rpkm $2.uniqs.recount -uniq -cache 1  $models $replacemodels"
python $ERANGEPATH/geneMrnaCountsWeighted.py $1 $2.rds $2.firstpass.rpkm $2.uniqs.recount -uniq -cache 1  $models $replacemodels

# count splice reads
echo "python $ERANGEPATH/geneMrnaCounts.py $1 $2.rds $2.splices.count -splices -noUniqs -cache 1  $models $replacemodels"
python $ERANGEPATH/geneMrnaCounts.py $1 $2.rds $2.splices.count -splices -noUniqs -cache 1  $models $replacemodels

# Alternative 1: find new regions outside of gene models with reads piled up 
echo "python $ERANGEPATH/findall.py RNAFAR $2.rds $2.newregions.txt -RNA -minimum 1 -nomulti -flag NM -log rna.log -cache 1"
python $ERANGEPATH/findall.py RNAFAR $2.rds $2.newregions.txt -RNA -minimum 1 -nomulti -flag NM -log rna.log -cache 1

# Alternative 1: filter out new regions that overlap repeats more than a certain fraction
echo "python $ERANGEPATH/checkrmask.py $3 $2.newregions.txt $2.newregions.repstatus $2.newregions.good -startField 1 -cache 1"
python $ERANGEPATH/checkrmask.py $3 $2.newregions.txt $2.newregions.repstatus $2.newregions.good -log rna.log -startField 1 -cache 1

# map all candidate regions that are within a given radius of a gene in bp
echo "python $ERANGEPATH/getallgenes.py $1 $2.newregions.good $2.candidates.txt -radius $4 -trackfar -cache  $models $replacemodels"
python $ERANGEPATH/getallgenes.py $1 $2.newregions.good $2.candidates.txt -radius $4 -trackfar -cache  $models $replacemodels

# make sure candidates.txt file exists
echo "touch $2.candidates.txt"
touch $2.candidates.txt

# calculate expanded exonic read density
echo "python $ERANGEPATH/normalizeExpandedExonic.py $1 $2.rds $2.uniqs.recount $2.splices.count $2.expanded.rpkm $2.candidates.txt $2.accepted.rpkm -cache  $models $replacemodels"
python $ERANGEPATH/normalizeExpandedExonic.py $1 $2.rds $2.uniqs.recount $2.splices.count $2.expanded.rpkm $2.candidates.txt $2.accepted.rpkm -cache  $models $replacemodels

# weigh multi-reads
echo "python $ERANGEPATH/geneMrnaCountsWeighted.py $1 $2.rds $2.expanded.rpkm $2.multi.count -accept $2.accepted.rpkm -multi -cache 1  $models $replacemodels"
python $ERANGEPATH/geneMrnaCountsWeighted.py $1 $2.rds $2.expanded.rpkm $2.multi.count -accept $2.accepted.rpkm -multi -cache 1  $models $replacemodels

# calculate final exonic read density
echo "python $ERANGEPATH/normalizeFinalExonic.py $2.rds $2.expanded.rpkm $2.multi.count $2.final.rpkm -multifraction -withGID -cache"
python $ERANGEPATH/normalizeFinalExonic.py $2.rds $2.expanded.rpkm $2.multi.count $2.final.rpkm -multifraction -withGID -cache

fi