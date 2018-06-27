#!/usr/bin/perl 
#use strict;
# read in HMM results and determines the start and end of enriched regions.

##### modified 04/08/08
##### 1. change allregions from .total.txt information to hmmout total consecutive regions.
##### 2. add minimum region length threshold, 50 bp  by default.
$minregionlength = 2; #2 * 25 = 50 bps
#$windowSize = 25;

$argc = @ARGV;
$argc == 6 || die "Provide probability threshold, significance level, window size, parameter filename, input filename and result file name.\n";

$threshold = $ARGV[0];
$siglevel = $ARGV[1];
$windowSize = $ARGV[2];
$parasfilename = $ARGV[3];
$inputfilename = $ARGV[4];
$outputfilename = $ARGV[5];

open(PARASFILE, $parasfilename);          
$datline = <PARASFILE>;
chop($datline);
@datas = split(" ", $datline);
$fgread = $datas[0];
for($j=0;$j<4;$j++)
{
        $datline = <PARASFILE>;
}
chop($datline);
@datas = split(" ", $datline);
$allregions = $datas[0];
$datline = <PARASFILE>;
chop($datline);
@datas = split(" ", $datline);
$fragmentWidth = $datas[0];

close(PARASFILE);
$lambda = $fgread*$fragmentWidth/(3100*0.9);           

$total = 0;
for($j=1;$j<25;$j++)
{
	open(INPUTFILE,$inputfilename.".chr".$j.".txt");
	my $count;
	my $first = 1;
	while($datline = <INPUTFILE>)
	{
		chop($datline);
	        @datas = split(" ", $datline);
	        $order = $datas[0];
		$proba = $datas[2];
		if($proba >= $threshold)
		{
 			if($first == 1)
			{
				$start = $order;
				$end = $order;
				$first = 0;
				$length = 1;
				$count = 1;
			}
			else
			{
				if($order == ($end +1))
				{# continue
					$end ++;
					$length ++;
				}
				else
				{# start a new region
					$start = $order;
					$end = $order;
        	                        $length = 1;
                	                $count ++;
				}
			}
		}
	}
	$total = $total + $count;
}
close(INPUTFILE);
print "there are $total regions.\n";
$allregions = $total;
$probcutoff = $siglevel/$allregions;   
$current = 0;
$prob = exp(-$lambda);
while($prob< (1-$probcutoff))
{
	$current ++;
	$logsum = 0;
        for($m=1;$m<=$current;$m++)
       	{
		$logsum = $logsum + log($lambda)-log($m);
	}
	$prob = $prob + exp($logsum-$lambda);
#	print "current = $current\tproba = $prob.\n";
}

$heightthreshold = $current;    
print "lambda = $lambda\theight threshold = $heightthreshold.\n";

open(OUTPUTFILE, ">$outputfilename.allregions.txt");
$total = 0;
for($j=1;$j<25;$j++)
{
	open(INPUTFILE,$inputfilename.".chr".$j.".txt");
	my $count;
	my $first = 1;
	while($datline = <INPUTFILE>)
	{
		chop($datline);
	        @datas = split(" ", $datline);
	        $order = $datas[0];
		$proba = $datas[2];
		$height = $datas[3];
		if($proba >= $threshold)
		{
 			if($first == 1)
			{
				$start = $order;
				$end = $order;
				$maxheight = $height;
				$first = 0;
				$length = 1;
				$count = 1;
			}
			else
			{
				if($order == ($end +1))
				{# continue
					$end ++;
					$length ++;
					if($height >$maxheight)
					{
						$maxheight = $height;
					}
				}
				else
				{# start a new region
					if(($maxheight > $heightthreshold)&&($length > ($minregionlength - 1)))
					{
                                                print OUTPUTFILE "$j\t",$start*$windowSize+1,"\t",($end+1)*$windowSize,"\t",$length*$windowSize,"\t$maxheight\n";
					}
					$start = $order;
					$end = $order;
					$maxheight = $height;
        	                        $length = 1;
                	                $count ++;
				}
			}
		}
	}
	if($maxheight> $heightthreshold)
        {	
                print OUTPUTFILE "$j\t",$start*25+1,"\t",($end+1)*25,"\t",$length*25,"\t$maxheight\n";
	}
	$total = $total + $count;
}
close(INPUTFILE);
close(OUTPUTFILE);
print "There are $total regions.\n";

sub decide()
{
	local($mylambda,$myprobcutoff)= @_;
	my @proba;
	my $j;
	my $k;
	my $m;
	my $middle = 30;
        my $logsum;
	$proba[0] = exp(-$mylambda);
	for($j = 1;$j<$middle;$j++) 	               
        {
		$logsum =0;
		for($m=1;$m<=$j;$m++)                
        	{                
        		$logsum = $logsum + log($lambda)-log($m);                
        	}
		$proba[$j] = exp(-$lambda+$logsum);
	}  
	for($j = $middle;$j<100;$j++)
        {
		$proba[$j] = 0;
	}
	$sum =0;
	for($j=0;$j<$middle;$j++)
	{
		$sum = $sum + $proba[$j];
	}
	my $sum = 0;   
	for($m=0;$m<=$middle;$m++)
        {
       		$sum = $sum + $proba[$m]*$proba[$m];                
        }
	my $trace = $sum;
	my $result = -1;
	if($trace>(1 - $myprobcutoff))
	{
		$result = 1;
	}
	else
	{
		$sum=0;
        	for($k=1;$k<=$middle;$k++)# just a large number to garantee achieve the significance level.
        	{
                	for($m=0;$m<=$middle;$m++)
                	{
                        	$sum = $sum + $proba[$m]*$proba[$m+$k];
                	}
			if((2*$sum+$trace) > (1-$myprobcutoff))
			{	
				$result = $k;
				last;
			}
		} 
        }
	$result;
}
