# an example shell script to combine multiple region calls into one partition
#

if [ -z "$1" ]; then
    PARTNAME=comb
else
	PARTNAME=$1
fi

if [ -z "$2" ]; then 
	MINSIZE=400
else
	MINSIZE=$2
fi

N=0
if [ $# -lt 2 ]; then
	echo 'usage: partition.sh name minSize datalist.file'
	echo
	echo 'where the datalist file is a list of region files'
	echo
else
	while read line
	do
		if [ $N -lt 1 ]; then
			FILELIST=''
		else
			FILELIST=$FILELIST,
		fi
		FILELIST=$FILELIST$line
		let N=N+1
	done < $3
	python $ERANGEPATH/partition.py $N.way $FILELIST $PARTNAME$N.part -minFeature $MINSIZE -nomerge -locid -norandom
fi
