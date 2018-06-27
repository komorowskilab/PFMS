#!/usr/bin/perl -w
#use strict;
#Written by: Steve Qin and Jianjun Yu 
# 05/20/08

use Cwd qw(realpath cwd);
use File::Basename;
use Getopt::Long;
Getopt::Long::Configure ('bundling_override');
$|=1;
my($fn, $dir) = fileparse(realpath($0));
my $pwd=cwd();
my $cmd=$0." ".join(" ",@ARGV); ####command line copy
my ($format,$help,$inputfilename,$inputcontrolfilename,$fmin,$fmax, $windowSize, $siglevel, $wigswitch,$seqswitch, $annotation, $name);

GetOptions('h'=>\$help,'format:s'=>\$format,'t:s'=>\$inputfilename,'c:s'=>\$inputcontrolfilename,'fmin:i'=>\$fmin,'fmax:i'=>\$fmax,'w|window:i'=>\$windowSize,
			's|sig:f'=>\$siglevel,'wig'=>\$wigswitch,'seq'=>\$seqswitch,'ann'=>\$annotation,'n|name:s'=>\$name);

if(defined $help || !$format || !$inputfilename) {
	print <<DES;
Usage: $fn <-format FORMAT -t TFILE -n NAME> [Options]

Options:
  -h               Show this help message
  -format          Format of tag file, "ELAND", "BED" or "Custom". REQUIRED
                   When choosing "Custom", the following info needs to be 
                   specified in the order:
                   1. The column to store genome file where match was found
                   2. The column to store base position of match
                   3. The column to store direction of match
                   The column numbers start at 1 and seperate with comma.
                   Example: the ELAND format can be presented as -format custom[7,8,9]
  -t TFILE         Treatment file name. REQUIRED
  -n, -name        Experiment name used to generate output file. REQUIRED
  -c CFILE         Control file name.
  -fmin            Minimal DNA fragment size. DEFAULT: 100
  -fmax            Maximal DNA fragment size. DEFAULT: 300
  -w, -window      Window size (bp). DEFAULT: 25
  -s, -sig         P-value threshold for peak detection. DEFAULT: 1e-3
  -wig             Whether to generate WIG file for UCSC genome brower
  -seq         	   Whether to extract peak sequences
  -ann         	   Whether to extract nearest gene information for peaks

Example: perl $fn -format ELAND -t stimu.inp -n stimu -fmin 100 -fmax 300 
	 perl $fn -format CUSTOM[6,7,8] -t chip.inp -c mock.inp -n chip-mock -fmin 100 -fmax 300 -w 25 -s 1e-3 -wig -seq -ann
	 
DES
exit;
}

$format=lc($format);
$fmin = 100 if(!$fmin);
$fmax = 300 if(!$fmax);
$windowSize = 25 if(!$windowSize);
$siglevel = 0.001 if(!$siglevel);
$fmin=int($fmin/$windowSize+0.5)*$windowSize;
$fmax=int($fmax/$windowSize+0.5)*$windowSize;
my $fragmentWidth = ($fmin+$fmax)/2;

open(SUMFILE, ">$name.sum");
my $logname = $name.".log";
open(LOGFILE, ">$logname");
print SUMFILE "program started at: ";    
$timestamp = localtime();
print "$timestamp\n"; 
print SUMFILE "$timestamp\n\n";

print SUMFILE "The command used is:\n";
print SUMFILE "$cmd\n\n";
print SUMFILE "Parameters specified:\n";
print SUMFILE "input treated name: $inputfilename\n";
print SUMFILE "treated sample data is in:\n";
open(INNAMEFILE, $inputfilename);
while($datline = <INNAMEFILE>)
{
        chop($datline);    
        if(length($datline) ==0)       
        {            
                next;           
        }
        print SUMFILE "$datline\n";       
}
close(INNAMEFILE);

if($inputcontrolfilename)
{
	print SUMFILE "input control name: $inputcontrolfilename\n";
	print SUMFILE "control sample data is in:\n";
	open(INNAMEFILE, $inputcontrolfilename); 
	while($datline = <INNAMEFILE>)
	{
        	chop($datline);          
        	if(length($datline) ==0)
        	{          
        	        next;   
        	}              
        	print SUMFILE "$datline\n";
	}
	close(INNAMEFILE);
}

print SUMFILE "output file name: $name\n";
print SUMFILE "format: $format\n";
print SUMFILE "fragement width: $fmin - $fmax bp\n";
print SUMFILE "window size: $windowSize\n";
print SUMFILE "significance level: $siglevel\n\n";

my $probcutoff = 1;##### can potentially be specified by the users.

my
 $tname=$cname=$name;
if($inputcontrolfilename) {
	$tname=$name."_treated";
	$cname=$name."_control";
}

##### first step, process raw data in treated sample and get summary.
print "##### process raw data in treated sample and get summary.\n";
print LOGFILE `perl $dir/summary.pl $inputfilename $format $fragmentWidth $siglevel $tname`;
die "Application Error: child process (summary.pl) exited abnormally.\n" if($?>>8 >0);
print SUMFILE "in treated sample:\n";
open(INFOFILE,"$tname.total.txt");
while($datline = <INFOFILE>)
{
	print SUMFILE "$datline";
}
print SUMFILE "\n";
close(INFOFILE);
print "done.\n";
$timestamp = localtime();
print "$timestamp\n";

##### second step, get genome-wide coverage profile.
print "##### get genome-wide coverage profile.\n";
print LOGFILE `perl $dir/chrwisewindow.pl $inputfilename $format $fmin $fmax $windowSize $tname`;
die "Application Error: child process (chrwisewindow.pl) exited abnormally.\n" if($?>>8 >0);
print "done.\n";
$timestamp = localtime();
print "$timestamp\n";

##### third step, use HMM to get enriched regions.
if(!$inputcontrolfilename)
{
	print "##### use HMM to get enriched regions.\n";
 	print LOGFILE `perl $dir/run.pl $tname $windowSize $siglevel`;
	die "Application Error: child process (run.pl) exited abnormally.\n" if($?>>8 >0);
	die "Child process (run.pl) died.\n" if($?>>8 >0);
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
}
##### fourth step, if there is control, get information.
if($inputcontrolfilename)
{
	print "##### process raw data in control sample and get summary.\n";
	print LOGFILE `perl $dir/summary.pl $inputcontrolfilename $format $fragmentWidth $siglevel $cname`;
	die "Application Error: child process (summary.pl) exited abnormally.\n" if($?>>8 >0);
	print SUMFILE "in control sample:\n";
	open(INFOFILE,"$cname.total.txt");
	while($datline = <INFOFILE>)
    	{
		print SUMFILE "$datline";
	}
	print SUMFILE "\n";
	close(INFOFILE);
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print "##### get genome-wide coverage profile.\n";
	print LOGFILE `perl $dir/chrwisewindow.pl $inputcontrolfilename $format $fmin $fmax $windowSize $cname`;
	die "Application Error: child process (chrwisewindow.pl) exited abnormally.\n" if($?>>8 >0);
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print "##### get summary of differences between treated and control samples.\n"; 
	print LOGFILE `perl $dir/summary.minus.pl $inputfilename $inputcontrolfilename $format $fragmentWidth $name`;
	die "Application Error: child process (summary.minus.pl) exited abnormally.\n" if($?>>8 >0);
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
   	print "##### get genome-wide coverage profile for the differences.\n";
	print LOGFILE `perl $dir/minus.pl $cname $tname $probcutoff $name`;
	die "Application Error: child process (minus.pl) exited abnormally.\n" if($?>>8 >0);
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print "##### use HMM to get enriched regions.\n";
	print LOGFILE `perl $dir/minus.run.pl $name $windowSize $siglevel`;
	die "Application Error: child process (minus.run.pl) exited abnormally.\n" if($?>>8 >0);
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
}
##### fifth step, find out number of region and total coverage.
open(REGIONFILE, "$name.allregions.txt");
my $count =0;
my $totalsize = $minwidth = $maxwidth = $minheigh = $maxheigh = 0;
while($datline = <REGIONFILE>)      
{
    chop($datline);
    my @datas = split(" ", $datline);
    my $start = $datas[1];
    my $end = $datas[2];
    my $width = $datas[3];
    my $heigh = $datas[4];
    if(($count==0)||($width > $maxwidth))
    {
        $maxwidth = $width;    
    }
    if(($count==0)||($width < $minwidth))
    {
        $minwidth = $width;
    }
    if(($count==0)||($heigh > $maxheigh))
    {
        $maxheigh = $heigh;
    }
    if(($count==0)||($heigh < $minheigh))
    {
        $minheigh = $heigh;
    }
    $count ++;
    $totalsize = $totalsize + $end-$start+1;
}

close(REGIONFILE);
$totalsize = $totalsize/1000000;
print SUMFILE "Number of enriched regions is $count.\n";
print SUMFILE "Total length of genome covered by the enriched regions: $totalsize Mb.\n";
print SUMFILE "size of the enriched regions range from $minwidth to $maxwidth.\n";
print SUMFILE "coverage for the enriched regions ranges from $minheigh to $maxheigh.\n";  
print SUMFILE "\nResult files:\n";
print SUMFILE "enriched regions are stored in $name.allregions.txt\n";

if(defined $wigswitch)
{
	if(!$inputcontrolfilename) {
        print "##### generate wig file for enriched regions.\n";
        print LOGFILE `perl $dir/wig.region.pl $name.allregions.txt $name $name $windowSize`;
		die "Application Error: child process (wig.region.pl) exited abnormally.\n" if($?>>8 >0);

	} else {
        print "##### generate wig file for enriched regions between treated and control samples.\n";
        print LOGFILE `perl $dir/wig.minus.region.pl $name.allregions.txt $tname $cname $name $windowSize`;
		die "Application Error: child process (wig.minus.region.pl) exited abnormally.\n" if($?>>8 >0);

	}
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print SUMFILE "WIG file of enriched regions is in $name.wig\n";
}

if(defined $seqswitch)
{
	print "##### obtain DNA sequences for enriched regions.\n";		
	print LOGFILE `perl $dir/get_seq.pl $name.allregions.txt`;
	die "Application Error: child process (get_seq.pl) exited abnormally.\n" if($?>>8 >0);
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print SUMFILE "FASTA format sequences in enriched regions are stored in $name.seq\n";         
}

if(defined $annotation)
{
	print "##### produce detailed information for enriched regions.\n";
	print "perl $dir/peakann.pl -a -o $name.annotation.txt $name.allregions.txt\n";
	print LOGFILE `perl $dir/peakann.pl -a -o $name.annotation.txt $name.allregions.txt`;
	die "Application Error: child process (peakann.pl) exited abnormally.\n" if($?>>8 >0);
	print "done.\n";
	$timestamp = localtime();
	print "$timestamp\n";
	print SUMFILE "Detailed annotation of enriched regions is in $name.annotation.txt\n";         
}
close(LOGFILE);

#### sixth step, clean up intermediate files, output wig, seq and annotation file if required.
print "##### clean up intermediate files.\n";
unlink <$tname.windowhitscount.chr*.txt>;
unlink <$tname.select.chr*.txt>;
unlink <$tname.location.txt>;
unlink <$tname.paras.txt>;
unlink <$tname.total.txt>;
unlink <$cname.windowhitscount.chr*.txt> if($inputcontrolfilename);
unlink <$cname.select.chr*.txt> if($inputcontrolfilename);
unlink <$cname.location.txt> if($inputcontrolfilename);
unlink <$cname.paras.txt> if($inputcontrolfilename);
unlink <$cname.total.txt> if($inputcontrolfilename);   
unlink <$name.paras.txt> if($inputcontrolfilename);
unlink <$name.hmmout.chr*.txt>; 
unlink <$name.dif.chr*.txt>; 

print "done.\n";
$timestamp = localtime();
print "$timestamp\n";
print SUMFILE "\nprogram ended at: ";
$timestamp = localtime(); 
print "$timestamp\n";
print SUMFILE "$timestamp\n";
close(SUMFILE);

