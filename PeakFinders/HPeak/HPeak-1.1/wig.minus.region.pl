#!/usr/bin/perl
use POSIX;
# get WIG file for enriched regions.
# Steve Qin 03/25/08

$argc = @ARGV;
$argc == 5 || die "Provide input region, input window hits count file names (case, control), output WIG file name and the window size 
(default is 25).\n";
$inputregionname = $ARGV[0];
$inputwindowhitscountname = $ARGV[1];
$controlwindowhitscountname = $ARGV[2];
$outputwigname = $ARGV[3];
$windowsize = $ARGV[4];
open(INPUTREGIONFILE, $inputregionname);
open(CONTROLWINDOWHITSCOUNTFILE, $controlwindowhitscountname);
open(OUTPUTWIGFILE, ">".$outputwigname.".wig");

$count = 0;
my @chr;
my @start;
my @width;
while($datline = <INPUTREGIONFILE>)
{
        chop($datline);
        @items = split(" ",$datline);
        $chr[$count] = $items[0];
        $start[$count]  = $items[1];
	$width[$count] = ($items[2] - $items[1] + 1)/25;
	$count++;
}
close(INPUTREGIONFILE);
print "There are $count regions.\n";
$allregions = $count;

##### output to WIG file
print OUTPUTWIGFILE "track name=\"$inputwindowhitscountname\" description=\"$inputwindowhitscountname\" type=\"wiggle_0\" color=50,50,150 yLineMark=0.0 yLineOnOff=on visibility=2\n";

$current = 0;
$inpeak = 0;
$counter = 0;
for($j=1;$j<25;$j++)
{
	open(INPUTWINDOWHITSCOUNTFILE, $inputwindowhitscountname.".windowhitscount.chr$j.txt");
        open(CONTROLWINDOWHITSCOUNTFILE, $controlwindowhitscountname.".windowhitscount.chr$j.txt");
	$chrnum = $j;
	if($j == 23)
	{
		$chrchr = "X";
	}
	elsif($j ==24)
	{
		$chrchr = "Y"; 
	}
	else
	{
		$chrchr = $j;
	}
	print "Chromosome $j.\n";
	while($datline = <INPUTWINDOWHITSCOUNTFILE>)
	{
		if($current == $allregions)
		{
			last;
		}
        	chop($datline);
        	@items = split(" ",$datline);
        	$chrstart  = $items[1];
		$peakheight = $items[3];
		$datline = <CONTROLWINDOWHITSCOUNTFILE>;
                chop($datline);
                @items = split(" ",$datline);
                $matchcontrolheight = $items[3];
		if($inpeak == 0)
		{
			if(($chr[$current] == $chrnum)&&($start[$current] == $chrstart))
        	        {
				print OUTPUTWIGFILE "fixedStep chrom=chr$chrchr start=$chrstart step=$windowsize\n";
				$inpeak = 1;
				print OUTPUTWIGFILE $peakheight - $matchcontrolheight,"\n";
				$counter = 1;
			}
		}
		else
		{	
			if($counter < $width[$current])
			{
                        	print OUTPUTWIGFILE $peakheight - $matchcontrolheight,"\n";
				$counter ++;      
			}
			else
			{
				$inpeak = 0;
				$current ++;
			}
                }		
	}
	close(INPUTWINDOWHITSCOUNTFILE);
        close(CONTROLWINDOWHITSCOUNTFILE);
}
close(OUTPUTWIGFILE);
