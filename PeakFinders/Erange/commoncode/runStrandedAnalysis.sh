#!/bin/bash
#
# runStrandedAnalysis.sh
# ENRAGE
#
# example: . ../commoncode/runStrandedAnalysis.sh mouse c2c12rna ../mm9repeats/rmask.db 20000
#       
#          assuming that we have rds database with the prefix c2c12rna.24R. 

# set ERANGEPATH to the absolute or relative path to ERANGE, if it's not in the environment

if [ -z "$ERANGEPATH" ]
then
    ERANGEPATH='../commoncode'
fi

echo 'runStrandedAnalysis.sh: version 4.1'

if [ -z "$1" ]
then
    echo
    echo 'usage:runStrandedAnalysis.sh genome rdsprefix repeatmaskdb bpradius'
    echo
    echo 'where rdsprefix is the name of the rds file without the .rds extension'
    echo 'use "none" for the repeatmaskdb if you do not have one'
    echo
else

# log the parameters
arguments=$1' '$2' '$3' '$4
echo 'running with settings: ' $arguments
python $ERANGEPATH/recordLog.py rna.log runStrandedAnalysis.sh "with parameters: $arguments"

# count the unique reads falling on the gene models ; the nomatch files are 
# mappable reads that fell outside of the Cistematic gene models and not the 
# unmappable of Eland (i.e, the "NM" reads)
python $ERANGEPATH/geneMrnaCounts.py $1 $2.rds $2.uniqs.count -stranded -markGID -cache 1

# calculate a first-pass RPKM to re-weigh the unique reads,
# using 'none' for the splice count
python $ERANGEPATH/normalizeExpandedExonic.py $1 $2.rds $2.uniqs.count none $2.firstpass.rpkm -cache

# recount the unique reads with weights calculated during the first pass
python $ERANGEPATH/geneMrnaCountsWeighted.py $1 $2.rds $2.firstpass.rpkm $2.uniqs.recount -stranded -uniq -cache 1

# count splice reads
python $ERANGEPATH/geneMrnaCounts.py $1 $2.rds $2.splices.count -stranded -splices -noUniqs -cache 1

# find new regions outside of gene models with reads piled up 
python $ERANGEPATH/findall.py RNAFARP $2.rds $2.newregions.txt -RNA -minimum 1 -nomulti -flag NM -strandfilter plus -log rna.log -cache 1
python $ERANGEPATH/findall.py RNAFARM $2.rds $2.newregions.txt -RNA -minimum 1 -nomulti -flag NM -strandfilter minus -log rna.log -cache 1 -append

# filter out new regions that overlap repeats more than a certain fraction
python $ERANGEPATH/checkrmask.py $3 $2.newregions.txt $2.newregions.repstatus $2.newregions.good -startField 1 -log rna.log -cache 1

# Alternative 2: use a precomputed list of "new" regions (outside of gene models)
#python $ERANGEPATH/regionCounts.py $3 $2.nomatch.bed $2.newregions.good $2.stillnomatch.bed
#python $ERANGEPATH/regionCounts.py $3 $2.rds $2.newregions.good 

# map all candidate regions that are within a given radius of a gene in bp
python $ERANGEPATH/getallgenes.py $1 $2.newregions.good $2.candidates.txt -radius $4 -trackfar -stranded -cache

# calculate expanded exonic read density
python $ERANGEPATH/normalizeExpandedExonic.py $1 $2.rds $2.uniqs.recount $2.splices.count $2.expanded.rpkm $2.candidates.txt $2.accepted.rpkm -cache

# weigh multi-reads
python $ERANGEPATH/geneMrnaCountsWeighted.py $1 $2.rds $2.expanded.rpkm $2.multi.count -accept $2.accepted.rpkm -stranded -multi -cache 1

# calculate final exonic read density
python $ERANGEPATH/normalizeFinalExonic.py $2.rds $2.expanded.rpkm $2.multi.count $2.final.rpkm -multifraction -withGID -cache

fi
