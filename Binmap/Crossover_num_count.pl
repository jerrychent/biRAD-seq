#!/usr/bin/perl -w
use strict;
#usage: perl $0 -i Offspring_SEG-map_BinFile.list -n 12 -o output
#usage: perl $0 -i /path/to/Offspring_SEG-map_BinFile.list -n 12 -o /path/to/output
#-n 表示染色体数目
use File::Basename;
use List::Util qw/sum/;
use Getopt::Std;
use vars qw($opt_i $opt_n $opt_o);
getopts ('i:n:o:');
my($chr,@C,$chro,$file,@list,%hash,$chromo,@number,$numb,$total);
open OUT,">$opt_o";
print OUT "#==========统计每个后代单株、每条染色体上重组断点的数目以及此单株总的重组断点数目==========\n";
for(1..$opt_n){
	$chr="Chr".sprintf("%02d",$_);
	push(@C,$chr);
}
$chro=join("\t",@C);
print OUT "#Sample_ID\t$chro\tSum\n";
open IN,$opt_i;##116-trans.rlt.bin.list
while(<IN>){chomp;
	$file=basename $_;
	open IN1,"$_";
	while(<IN1>){chomp;
		@list=split;
		$hash{$list[0]}++;
	}
	close IN1;		
	for(1..$opt_n){
		$chromo="chromosome".sprintf("%02d",$_);
		push(@number,$hash{$chromo}-1);
	}
	$numb=join("\t",@number);
	$total=(sum @number);
	print OUT "$file\t$numb\t$total\n";
	@number=();
	%hash=();
}
close IN;
close OUT;
