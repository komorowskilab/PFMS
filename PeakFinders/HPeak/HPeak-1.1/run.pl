#!/usr/bin/perl
#use strict;
# Steve Qin 11/30/07
use Cwd qw(realpath cwd);
use File::Basename;
my($fn, $dir) = fileparse(realpath($0));

$argc = @ARGV;
$argc == 3 || die "Provide input file name, window size and significance level.\n";

$name = $ARGV[0];
$windowSize = $ARGV[1];
$siglevel = $ARGV[2];
open(SIZEINFOFILE,"$dir/data/chromsizes.txt");
$count = 0;
my @sizes;
while($dataline = <SIZEINFOFILE>)
{
	chop($dataline);
	@datas = split(" ",$dataline);
	$sizes[$count] = $datas[0];
	$count ++;
}
for($k=1; $k<25; $k++)
{
	$com = "$dir/chiphmm ".$name.".select.chr".$k.".txt ".$name.".paras.txt ".$name.".hmmout.chr".$k.".txt $windowSize ".$sizes[$k-1];
	print "$com\n";
	system($com);
}
$com = "perl $dir/autoregion.pl 1 $siglevel $windowSize ".$name.".paras.txt ".$name.".hmmout ".$name;
print "$com\n";
system($com);

