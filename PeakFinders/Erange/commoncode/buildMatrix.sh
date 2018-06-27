#!/bin/bash
echo 'buildMatrix.sh: version 1.1'

indexPrev=0
indexCur=0

truncateRPKM=""
if [ $# -eq 3 ]; then
    truncateRPKM="-truncate "$3
fi

if [ $# -eq 4 ]; then
    truncateRPKM="-rescale -truncate "$3
fi

if [ $# -lt 2 ]; then
	echo 'usage: buildMatrix.sh name datalist.file [truncateRPKM] [-rescale]'
	echo
	echo 'where the datalist file is a comma-delimited list of prefix and rds-files'
	echo
else
	python $ERANGEPATH/recordLog.py buildMatrix.log buildMatrix.sh "with parameters: $1 $2 $truncateRPKM"
	while read line
	do 
		prefix=`echo $line | cut -f 1 -d ','`
		filename=$prefix.partcount
		if [ -e $filename ]; then
			if [ $indexCur -lt 1 ]; then
				echo "building $1.step0"
				echo -e '\t' > $1.step0
				cut -f 1 $filename >> $1.step0
				indexCur=1
			fi
				python $ERANGEPATH/buildMatrix.py $1.step$indexPrev $filename $1.step$indexCur $truncateRPKM
				rm $1.step$indexPrev
				let indexPrev=indexPrev+1
				let indexCur=indexCur+1
		else
			echo "could not find $filename - skipping"
			python $ERANGEPATH/recordLog.py buildMatrix.log buildMatrix.sh "could not find $rds - skipping"
		fi
	done < $2      
	mv $1.step$indexPrev $1.matrix.tab
fi
