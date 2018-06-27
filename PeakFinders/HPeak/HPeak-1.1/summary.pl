#!/usr/bin/perl
use POSIX;
use Cwd qw(realpath);
use File::Basename; 
# Steve Qin 08/01/08
$fragmentWidth = 300;

my $script_path=realpath($0);
my($ef, $dir) = fileparse($script_path);

$argc = @ARGV;
$argc == 5 || die "Provide input files names that include multiple realign.txt filenames, format, fragment size, significant level and result file name.\n";

$innamename = $ARGV[0];
open(INNAMEFILE, $innamename) || die "Can not open $innamename\n";
$format = $ARGV[1];
$format=lc($format);
$fragmentWidth = $ARGV[2];
$siglevel = $ARGV[3];
open(OUTFILE1, ">".$ARGV[4].".location.txt") || die "$!\n";
open(OUTSUMFILE,">".$ARGV[4].".total.txt") || die "$!\n";
open(OUTPARFILE,">".$ARGV[4].".paras.txt") || die "$!\n";

$count = 0;
my @filenames;
while($datline = <INNAMEFILE>)
{
    chomp($datline);
	if(length($datline) ==0)
	{
		next;
	}
	$filenames[$count] = $datline;
	$count ++;
}
$filecount = $count;
print "There are $filecount files.\n";

for($i = 0;$i < 24; $i++)
{
        $hitscount[$i] = 0;
	$posreadscount[$i] = 0;
	$negreadscount[$i] = 0;
}
$totalcount = 0;
$count = 0;
for($i=0;$i<$filecount;$i++)
{
	$currentfilename = $filenames[$i];
	&realignfileread($currentfilename,$format,*totalcount, *count,*hits, *hitscount,*posreadscount, *negreadscount);##### newly added 04/13/08
	print "total count = $totalcount.\n";
        print "count = $count.\n";
}
$totalread = $count;

$count = 0;
for($i = 0;$i < 24; $i++)
{
	if($hitscount[$i]==0)
        {
                next;
        }
	elsif($hitscount[$i]==1)
        {
                $count ++;
		next;
        }
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
	$count ++;
}
$allregion = $count;
print "total region = $count.\n";

$totalpeaks = 0;
$peakcover = 0;
$totalcover = 0;
$totalinpeakreads = 0;
my @peakwidth; 
my @inpeakreadscount;
for($i = 0;$i < 24; $i++)
{
	if($hitscount[$i]==0)
        {
		next;
	}
	my @pos;
	for($j=0;$j<$hitscount[$i];$j++)
	{
		$pos[$j] = $hits[$i][$j];
	}
        @sortpos = sort {$a <=> $b} @pos;
	$count = 1;
	$start = $sortpos[0];
	print OUTFILE1 $i + 1,"\t$start\t";
	if($hitscount[$i]==1)
	{
		$end = $sortpos[0] + $fragmentWidth -1;                 
	        $regionwid = $end - $start +1;
        	$totalcover = $totalcover + $regionwid; 
		print OUTFILE1 "\n";
		next;
	}
        for($j=1;$j<$hitscount[$i];$j++) 
        { 
		$dis = $sortpos[$j] - $sortpos[$j-1];
		if(($sortpos[$j] - $sortpos[$j-1])< ($fragmentWidth))
		{#extension continues
	                print OUTFILE1 "$sortpos[$j]\t";
			$count ++;
			next;
		}
		else 
		{#extension stops
			print OUTFILE1 "\n";
			$end = $sortpos[$j-1] + $fragmentWidth -1;
			$regionwid = $end - $start +1;
			$totalcover = $totalcover + $regionwid;
			$lambda = $totalread /(0.9*3100000000) * $regionwid;
			$sum = 1;
			for($k=1;$k<$count;$k++)
			{
				$logsum = 0;
				for($m=1;$m<=$k;$m++)
                        	{
					$logsum = $logsum + log($lambda)-log($m);
				}
				$sum = $sum + exp($logsum);
			}
			$prob = 1 - exp(-$lambda)* $sum;
			if($prob < 0.0000000001)
			{
				$prob = 0;
			}
                        if($prob < $siglevel/$allregion)
			{
				$totalinpeakreads = $totalinpeakreads + $count;
				$peakwidth[$totalpeaks] = $regionwid;
                                $inpeakreadscount[$totalpeaks] = $count;
				$totalpeaks ++;
				$peakcover = $peakcover + $regionwid;
			}
        	        $start = $sortpos[$j];
		        print OUTFILE1 $i + 1,"\t$start\t";
                        $count = 1;
		}
	} 
        $end = $sortpos[$hitscount[$i]-1] + $fragmentWidth -1;
        $regionwid = $end - $start +1;
	$totalcover = $totalcover + $regionwid;
        $lambda = $totalread /(0.9*3100000000) * $regionwid;
        $sum = 1;
        for($k=1;$k<$count;$k++)
        {
                $logsum = 0;
                for($m=1;$m<=$k;$m++)
                {
                        $logsum = $logsum + log($lambda)-log($m);
                }
                $sum = $sum + exp($logsum);
        }
        $prob = 1 - exp(-$lambda)* $sum;
        if($prob < 0.0000000001)
        {
                $prob = 0;
        }
        if($prob < $siglevel/$allregion)
        {
                $totalinpeakreads = $totalinpeakreads + $count;
                $peakwidth[$totalpeaks] = $regionwid;
                $inpeakreadscount[$totalpeaks] = $count;
                $totalpeaks ++;
                $peakcover = $peakcover + $regionwid;
         }
   print OUTFILE1 "\n";
}
$avecover = $totalcover /$totalread;
$totalcover = $totalcover /1000000;
print "total coverage is $totalcover MB.\n";
print "average coverage per read is $avecover KB.\n";
close(OUTFILE1);

@sortpeakwidth = sort {$a <=> $b} @peakwidth;
$medianwidth = $sortpeakwidth[floor(0.5*$totalpeaks)];

print "total number of reads = $totalcount.\n";
print "total number of successfully aligned reads = $totalread. (",100*$totalread/$totalcount," %)\n";
print OUTSUMFILE "Total sequenced reads (millions)\t",$totalcount/1000000,"\n";
print OUTSUMFILE "Total uniquely mapped reads (millions)\t",$totalread/1000000," (",100*$totalread/$totalcount," %)\n";
close(OUTSUMFILE); 

print OUTPARFILE $totalread/1000000,"\n";
print OUTPARFILE $totalinpeakreads/1000000,"\n";
print OUTPARFILE $peakcover/1000000,"\n";
print OUTPARFILE "$medianwidth\n";
print OUTPARFILE "$totalpeaks\n";
print OUTPARFILE $fragmentWidth,"\n";#to accout for linear descreasing at the tail.
close(OUTPARFILE);

sub realignfileread {
         local ($subcurrentfilename, $subformat, *subtotalcount, *subcount, *subhits, *subhitscount, *subposreadscount, *subnegreadscount) =@_;
         if($subcurrentfilename =~ /gz$/)
         {
                 system("gunzip ".$subcurrentfilename);
                 $namelength = length($subcurrentfilename);
                 $newname = substr($subcurrentfilename,0,$namelength-3);
                 open(INREADFILE, $newname) || die "$!\n";
         }
         elsif($subcurrentfilename =~ /zip$/)
         {
                 system("unzip ".$subcurrentfilename);
                 $newname = substr($subcurrentfilename,0,$namelength-4);
                 open(INREADFILE, $newname.".txt") || die "$!\n";
         }
         else
         {
                 open(INREADFILE, $subcurrentfilename) || die "$!\n";
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
    } 
    elsif ($subformat=~/custom\[(.*)\]/)
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
         $subposreadscount[$chr-1] ++;
    }
    elsif($direct eq "-" || $direct eq "R")
    {
         $subhits[$chr-1][$subhitscount[$chr-1]] = $posit-$fragmentWidth +1;
         $subnegreadscount[$chr-1] ++;
    }
    $subhitscount[$chr-1] ++;
    $subcount++;
}
close(INREADFILE);
print "subcount = $subcount.\n";
}

