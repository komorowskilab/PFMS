#!/usr/bin/perl
use POSIX;
# Steve Qin 12/03/07

$argc = @ARGV;
$argc == 4 || die "Provide unsimulated, simulated data files, threshold and result file name.\n";

$unsimfilename = $ARGV[0];
$simulfilename = $ARGV[1];
$threshold = $ARGV[2];
$outputfilename = $ARGV[3];

open(PARASFILE,$outputfilename.".paras.txt");
$datline = <PARASFILE>;
$datline = <PARASFILE>;
$datline = <PARASFILE>;
$datline = <PARASFILE>;
chop($datline);             
@datas = split(" ", $datline);        
$fgrate = $datas[0];      
$datline = <PARASFILE>;
chop($datline);                 
@datas = split(" ", $datline);            
$bgrate = $datas[0];            
$folds = $fgrate/$bgrate;
print "$bgrate\t$fgrate\t$folds\n";               
close(PARASFILE);

for($j=1;$j<=24;$j++)
{
	print "Chromosome $j:\n";
	&subtract($j,$unsimfilename,$simulfilename,$outputfilename,$folds,$threshold);
}

sub subtract()
{
	local($mychr, $myunsimfilename,$mysimulfilename,$myoutputfilename, $myfolds, $mythreshold) = @_;
	$localunsimfilename = $myunsimfilename.".windowhitscount.chr".$mychr.".txt";
	$localsimulfilename = $mysimulfilename.".windowhitscount.chr".$mychr.".txt";
	$localoutputfilename = $myoutputfilename.".dif.chr".$mychr.".txt";
	open(UNSIMFILE, $localunsimfilename); 
	open(SIMFILE, $localsimulfilename); 
        open(OUTPUTFILE, ">".$localoutputfilename);
my $count = 0;
while($datline = <UNSIMFILE>)
{
        chop($datline);
	@datas = split(" ", $datline);
	
	$bghits = $datas[3];
	$datline = <SIMFILE>;
	chop($datline);
        @datas = split(" ", $datline);
        $hits = $datas[3];
	$dif = $hits - $bghits * $myfolds;	
	if($dif >= $mythreshold)
	{
		print OUTPUTFILE "$count\t$dif\t$hits\t",$bghits * $myfolds,"\n";
	}
	$count ++;
}
close(UNSIMFILE);
close(SIMFILE);
close(OUTPUTFILE);
}
