#!/usr/bin/perl
use POSIX;
# Steve Qin 11/02/07

$argc = @ARGV;
$argc == 5 || die "Provide input case and control files names that include multiple eland.result.txt filenames, input file format, fragment size, and result file name.\n";

$innamename = $ARGV[0];
open(INNAMEFILE, $innamename); 
$incontrolnamename = $ARGV[1];
open(INCONTROLNAMEFILE, $incontrolnamename);
$format = $ARGV[2];
$fragmentWidth = $ARGV[3];
open(OUTPARFILE,">".$ARGV[4].".paras.txt");

$count = 0;
my @filenames;
while($datline = <INNAMEFILE>)
{
        chop($datline);
        if(length($datline) ==0)
        {
                next;
        }
        $filenames[$count] = $datline;
	$count ++;
}
$filecount = $count;
print "There are $filecount case files.\n";
for($i = 0;$i < 24; $i++)
{
        $hitscount[$i] = 0;
}
$totalcount = 0;
$count = 0;
for($i=0;$i<$filecount;$i++)
{
	$currentfilename = $filenames[$i];
	&realignfileread($currentfilename,$format,*totalcount, *count,*hits, *hitscount);
	print "total count = $totalcount.\n";
        print "count = $count.\n";
}
$totalread = $count;
$count = 0;
my @controlfilenames;
while($datline = <INCONTROLNAMEFILE>)
{
        chop($datline);
        if(length($datline) ==0)
        {
                next;
        }
        $controlfilenames[$count] = $datline;
        $count ++;
}
$filecount = $count;     
print "There are $filecount control files.\n";
for($i = 0;$i < 24; $i++)
{
        $controlhitscount[$i] = 0;
}
$totalcontrolcount = 0;
$controlcount = 0;
for($i=0;$i<$filecount;$i++)
{
        $currentfilename = $controlfilenames[$i];
        &realignfileread($currentfilename,$format,*totalcontrolcount, *controlcount,*controlhits, *controlhitscount);
        print "total control count = $totalcontrolcount.\n";
        print "control count = $controlcount.\n";
}
$totalcontrolread = $controlcount; 

$ratio = $totalread/$totalcontrolread;
print "$totalread\t$totalcontrolread\t$ratio\n";

$count = 0;
for($i = 0;$i < 24; $i++)
{
        my @pos;
        for($j=0;$j<$hitscount[$i];$j++)
        {
                $pos[$j] = $hits[$i][$j];
        }
        @sortpos = sort {$a <=> $b} @pos;
        for($j=1;$j<$hitscount[$i];$j++)
        {
                if(($sortpos[$j] - $sortpos[$j-1])< ($fragmentWidth))
                {#extension continues
                        next;
                }
                else
                {#extension stops
                        $count ++;
                }
        }       
}
$allregions = $count;
print "total region = $count.\n";

$probcutoff = 0.001/$allregions;
print "$allregions\n";
$order = 1;
$totalpeaks = 0;
$peakcover = 0;
$totalinpeakreads = 0;
$totalinpeakcontrolreads = 0;
$totalinpeakdif = 0;
my @peakwidth; 
my @inpeakdifcount;
$badcount = 0;
for($i = 0;$i < 24; $i++)
{
	print "chr",$i+1,"\n";
	my @pos;
	for($j=0;$j<$hitscount[$i];$j++)
	{
		$pos[$j] = $hits[$i][$j];
	}
        @sortpos = sort {$a <=> $b} @pos;

        my @controlpos;
        for($j=0;$j<$controlhitscount[$i];$j++)
        {
                $controlpos[$j] = $controlhits[$i][$j];
        }
        @sortcontrolpos = sort {$a <=> $b} @controlpos;

	$count = 1;
	$start = $sortpos[0];
	my $j;
	$skipstart = 0;
        for($j=1;$j<$hitscount[$i];$j++) 
        { 
		if(($sortpos[$j] - $sortpos[$j-1])< ($fragmentWidth))
		{#extension continues
			$count ++;
			next;
		}
		else 
		{#extension stops
			$end = $sortpos[$j-1] + $fragmentWidth;
			$regionwid = $end - $start;

			$matchcontrolcount = &checkcontrol($start,$end,$controlhitscount[$i],*skipstart, *sortcontrolpos);
                        $baseline = $totalread /(0.9*3100000000);
			$dif = $count - $ratio * $matchcontrolcount;
			if($dif > 0)
			{
				if($dif > 100)
				{
					$outcome = 1;
				}
				else
				{
					$outcome = &decide($baseline, $regionwid, $dif, $probcutoff);
				}
	                        if($outcome ==1)
        	                {
					$totalinpeakreads = $totalinpeakreads + $count;
					$totalinpeakcontrolreads = $totalinpeakcontrolreads + $matchcontrolcount;
					$totalinpeakdif = $totalinpeakdif + $dif;
					$peakwidth[$totalpeaks] = $regionwid;
                        	        #$inpeakreadscount[$totalpeaks] = $count;
					$inpeakdifcount[$totalpeaks] = $dif;
					$totalpeaks ++;
					$peakcover = $peakcover + $regionwid;
				}
			}
        	        $start = $sortpos[$j];
			$order ++;
                        $count = 1;
		}
	} 
}
print "badcount = $badcount.\n";

@sortpeakwidth = sort {$a <=> $b} @peakwidth;
$medianPeakwidth = $sortpeakwidth[floor(0.5*$totalpeaks)];

print "total number of reads = $totalcount.\n";
print "total number of successfully aligned reads = $totalread. (",100*$totalread/$totalcount," %)\n";
print "total number of control reads = $totalcontrolcount.\n";
print "total number of successfully aligned control reads = $totalcontrolread. (",100*$totalcontrolread/$totalcontrolcount," %)\n";

print OUTPARFILE $peakcover/1000000,"\n";    
print OUTPARFILE $totalinpeakreads,"\n";
print OUTPARFILE $totalinpeakcontrolreads,"\n";    
print OUTPARFILE "$totalread\n";
print OUTPARFILE "$totalcontrolread\n";
print OUTPARFILE $fragmentWidth,"\n";# to account for linear descreasing at the tail.
print OUTPARFILE "$medianPeakwidth\n";
close(OUTPARFILE);

sub checkcontrol
{
	local($localstart, $localend,$localchrsize, *localskipstart,*localsortcontrolpos) = @_;	
	my $j;
	my $count = 0;
	my $newstart = $localskipstart;
	for($j=$localskipstart;$j<$localchrsize;$j++)
	{
		if(($localstart - $localsortcontrolpos[$j]) > $fragmentWidth)
		{
			$newstart = $j;
			next;
		}
		if($localend < $localsortcontrolpos[$j])
                {
			last;
		}
		if((($localstart - $localsortcontrolpos[$j]) <= $fragmentWidth)&&($localend >= $localsortcontrolpos[$j]))
		{ 
                        $count ++;
                }
	}
	$localskipstart = $newstart;
	$count;
}

sub realignfileread
{
	local ($subcurrentfilename, $subformat, *subtotalcount, *subcount, *subhits, *subhitscount) = @_; 
	if($subcurrentfilename =~ /gz$/)
	{
		system("gunzip ".$subcurrentfilename);
		$namelength = length($subcurrentfilename);		
		$newname = substr($subcurrentfilename,0,$namelength-3); 
                open(INREADFILE, $newname);
	}
	elsif($subcurrentfilename =~ /zip$/)
        {
                system("unzip ".$subcurrentfilename);
		$newname = substr($subcurrentfilename,0,$namelength-4);
                open(INREADFILE, $newname.".txt");
	}
	else
	{
		print "$subcurrentfilename.\n";
		open(INREADFILE, $subcurrentfilename);		
	}
while($datline = <INREADFILE>)
{
    my $chr=$posit=$direct=$loc=0;

    chomp($datline);
    my @datas=split /[\s\t]+/, $datline;
    $subtotalcount ++;

    if($subformat eq "eland")
    {
        if(($datline=~ /chr/)&&($datline =~/\.fa/))
        {
                $posit = $datas[7];
                $direct = $datas[8];
                $loc=$datas[6];
        }
    } elsif ($subformat=~/custom\[(.*)\]/)
    {
	$tempo = $1;
        if(($datline=~ /chr/)&&($datline =~/\.fa/))
        {
                my @infos=split /\,\s*/,$tempo;
                $loc=$datas[$infos[0]-1];
                $posit=$datas[$infos[1]-1];
                $direct=$datas[$infos[2]-1];
        }
    } elsif ($subformat eq "bed") {
         $chr =$2 if($datas[0]=~/^(chr)?(\d+|X|Y)$/);
         $posit = $datas[1]+1;
         $direct = $datas[3];
    }
    $chr=$1 if($loc=~/chr(\d+|X|Y)\.fa/ || $loc=~/chromosome\.(\d+|X|Y)\.fa/);
    next if(!$chr);
    $chr=23 if($chr eq "X");
    $chr=24 if($chr eq "Y");
    if($direct eq "+" || $direct eq "F")
    {
         $subhits[$chr-1][$subhitscount[$chr-1]] = $posit;
    }
    elsif($direct eq "-" || $direct eq "R")
    {
         $subhits[$chr-1][$subhitscount[$chr-1]] = $posit-$fragmentWidth +1;
         ##### newly added 04/13/08 START
         ##### newly added 04/13/08 END
    }
    $subhitscount[$chr-1] ++;
    $subcount++;
}
close(INREADFILE);
print "subcount = $subcount.\n";
}

sub decide()
{
	local($mybaseline, $myregionwid, $mydif, $myprobcutoff)= @_;
        my $lambda = $mybaseline * $myregionwid;
	my @proba;
	my $j;
	my $k;
	my $m;
	my $middle = 30;
	my $upper = 150; # need $upper - $middle > 100.
        my $logsum;
	$proba[0] = exp(-$lambda);
	for($j = 1;$j<$middle;$j++) 	               
        {
		$logsum =0;
		for($m=1;$m<=$j;$m++)                
        	{                
        		$logsum = $logsum + log($lambda)-log($m);                
        	}
		$proba[$j] = exp(-$lambda+$logsum);
	}  
	for($j = $middle;$j<$upper;$j++)
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
	$sum=0;
	my $finish = 0;
	my $result = -1;
        for($k=1;$k<=$mydif;$k++)
        {
                for($m=0;$m<=$middle;$m++)
                { 
                       $sum = $sum + $proba[$m]*$proba[$m+$k];
                }
		if((2*$sum+$trace) > (1-$myprobcutoff))
		{
			$result = 1;
			$finish = 1;
			last;
		} 
        }
	if(($finish == 0)&&((2*$sum+$trace) <= (1-$myprobcutoff)))
        {
        	$result = 0;
        }
	$result;
}
