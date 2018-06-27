#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   bar2wig.pl
# 
# Description:
#   Convert BAR files to WIG (BedGraph) files
# 
# Usage:
#   perl bar2wig.pl <*.bar> <*.wig>
# 
# Author:
#   Xi Wang, wang-xi05@mails.thu.edu.cn
# 
# Date:
#   Thu Jan 27 00:19:01 CST 2011
#
########################################

use strict;
my $usage = "$0 <*.bar> <*.wig>\n";
my $infile = shift || die $usage;
my $outfile = shift || die $usage;
open(IN, $infile) || die "Can't open $infile for reading!\n";
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";

my ($chr, $s, $score);
my $chrN_inserted = 0;
while(<IN>)
{
	chomp;
	my @col = split;
	$s = $col[1] - 1;
	$score = $col[-1];
	if(!$chrN_inserted){
		$chr = $col[0];
		print OUT "track type=wiggle_0 name=\"SeqSite TFBS\" description=\"binding sites detected by SeqSite\"\nvariableStep chrom=$chr span=1\n";
		$chrN_inserted =1
	}
	print OUT "$s\t$score\n";
}
close IN;
close OUT;
