#!/usr/bin/perl
#use strict;
# read in HMM results and determines the start and end of enriched regions.

$minregionlength = 2; #2 * 25 = 50 bps

$argc = @ARGV;
$argc == 6 || die "Provide probability threshold, significance level, window size, input filename, parameter file and result file name.\n";

$threshold = $ARGV[0];
$siglevel = $ARGV[1];
$windowSize = $ARGV[2];
$inputfilename = $ARGV[3];
$parasfilename = $ARGV[4];
$outputfilename = $ARGV[5];
my @chr;
my @startpos;
my @endpos;
my @regionwid;
my @supheight;

$totalcount = 0;
for($j=1;$j<25;$j++)
{
	open(INPUTFILE,$inputfilename.".chr".$j.".txt");
	my $count;
	my $first = 1;
	my $reached = 0;
	while($datline = <INPUTFILE>)
	{
		chop($datline);
	        @datas = split(" ", $datline);
	        $order = $datas[0];
		$proba = $datas[2];
		$height = $datas[3];
		if($proba >= $threshold)
		{
			$reached = 1;
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
					$chr[$totalcount] = $j;
					$startpos[$totalcount] = ($start-1)*$windowSize+1;
					$endpos[$totalcount] = $end*$windowSize;
					$regionwid[$totalcount] = $length *$windowSize;
					$supheight[$totalcount] = $maxheight;					
					$start = $order;
					$end = $order;
					$maxheight = $height;
        	                        $length = 1;
                	                $count ++;
					$totalcount ++;
				}
			}
		}
	}
	if($reached ==1)
	{
		$chr[$totalcount] = $j;
        	$startpos[$totalcount] = ($start-1)*$windowSize+1;
        	$endpos[$totalcount] = $end*$windowSize;          
        	$regionwid[$totalcount] = $length *$windowSize;
        	$supheight[$totalcount] = $maxheight;
		$totalcount ++;
		print "totalcount = $totalcount\n";
	}
}
close(INPUTFILE);

open(PARASFILE, $parasfilename);
for($j=0;$j<4;$j++)
{
	$datline = <PARASFILE>;           
}
chop($datline);
@datas = split(" ", $datline); 
$fgread = $datas[0];
$datline = <PARASFILE>;  
$datline = <PARASFILE>;
chop($datline);             
@datas = split(" ", $datline);
$fragmentWidth = $datas[0];

print "foreground read = $fgread\n";
close(PARASFILE);
$lambda = $fgread/1000000*$fragmentWidth/(3100*0.9);
$probcutoff = $siglevel/$totalcount;
$heightcutoff = &decide($lambda,$probcutoff);
print "There are $totalcount regions.\n";
print "lambda= $lambda.\n";
print "max height = $heightcutoff\n";
open(OUTPUTFILE, ">".$outputfilename.".allregions.txt");
for($j=0;$j<$totalcount;$j++)
{
	if(($supheight[$j] > $heightcutoff)&&($regionwid[$j] > ($minregionlength - 1)*$windowSize))
        {        
               print OUTPUTFILE "$chr[$j]\t$startpos[$j]\t$endpos[$j]\t$regionwid[$j]\t$supheight[$j]\n";
        }			
}
close(OUTPUTFILE);

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
	print "result= $result\n";
	$result;
}
