#!/bin/bash
echo 'regionCounts.sh: version 1.0'

cachepages=""
if [ $# -eq 3 ]; then
    cachepages="-cache "$3
fi

if [ $# -lt 2 ]; then
	echo 'usage: regionCounts.sh partitionfile datalist.file [cachevalue]'
	echo
	echo 'where the datalist file is a comma-delimited list of prefix and rds-files'
	echo
else
	arguments=$1' '$2' '$cachepages
	python $ERANGEPATH/recordLog.py regionCounts.log regionCounts.sh "with parameters: $arguments"
	while read line
	do 
		prefix=`echo $line | cut -f 1 -d ','`
		rds=`echo $line | cut -f 2 -d ','`
		if [ -e $rds ]; then
			python $ERANGEPATH/regionCounts.py $1 $rds $prefix.partcount -force -nomerge -rpkm $cachepages
		else
			echo "could not find $rds - skipping"
			python $ERANGEPATH/recordLog.py regionCounts.log regionCounts.sh "could not find $rds - skipping"
		fi
	done < $2      
fi
