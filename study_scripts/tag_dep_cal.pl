#!/usr/bin/perl -w
#usage :perl $0 rice26-AD2-1.sort.bam.dep.tags.reads
#
my ($out,@list, @set,$num);
open IN,$ARGV[0];
$out=$ARGV[0].".dep_num.txt";
open OUT,">./$out";
while(<IN>){chomp;
	@list=split;
	@set=split (/,/,$list[2]);
	$num=@set;
	print OUT "$list[0]\t$list[1]\t$num\n";
}
close IN;
close OUT;
