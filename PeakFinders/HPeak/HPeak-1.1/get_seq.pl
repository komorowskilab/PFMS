#!/usr/bin/perl
use Cwd qw(realpath cwd);
use File::Basename;
my($fn, $dir) = fileparse(realpath($0));

my $file=$ARGV[0];
open(IN,"<$file") || die "can not open the file\n";
open(OUT2,">$file.seq") || die "can not open the output seq file\n";
$count = 0;
while(my $line=<IN>) {
	chomp($line);
	split(" ",$line);
	my ($chr,$start,$end,$len,$freq)=($_[0],$_[1],$_[2],$_[2]-$_[1]+1,$_[3]); #for stat1, $freq should be $_[1]
	if(!($chr =~ /chr/))
	{
		$chr=lc("chr".$chr);
	}
	if($chr eq "chr23") {
		$chr="chrX";
	}
	elsif($chr eq "chr24") {
		$chr="chrY";
	}
	my ($good,$seq)=getSeq($chr,$start,$end,$len);
	next if($good !=1);
	$count++;
	print OUT2 $seq;
}
print "count = $count\n";

sub getSeq {
	local($chr,$start,$end,$len)=@_;
	open(IN1,"<$dir/data/chromFa/$chr.fa") || die "can not open chromosome seq file\n";
#	print "<$dir/data/chromFa/$chr\.fa\n";
	my $offset=1+length($chr);
	$start=1 if($start<1);
	my $buffer="";
	seek(IN1, $start+$offset+int($start/50),0);
	read(IN1,$buffer,$end-$start+1+int($end/50)-int($start/50));
	close(IN1);
	my $masked_seq=()=$buffer=~/N/g;
	print "width is 0, $chr,$start,$end,$len.\n" if(length($buffer)==0);
	return (0,0) if(($masked_seq/length($buffer))>0.8);#####
	$buffer=~s/\s//g;
	while($buffer=~/^[ACGT]{0,15}N+/ || $buffer=~/N+[ACGT]{0,15}$/) {
		$buffer=~s/^[ACGT]{0,15}N+//;
		$buffer=~s/N+[ACGT]{0,15}$//;
	}
	if(length($buffer)>=50) {#####>100
	return (1,">$chr:$start-$end;$freq\n".$buffer."\n");
	}
	else {
	return (0,0);
	}
}
