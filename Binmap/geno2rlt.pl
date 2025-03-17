#/usr/bin/perl -w
use strict;
#usage: perl $0 B4-1.ge.geno  B4-1.ge.geno.rlt
my (@list);
open IN,$ARGV[0];
open OUT,">./$ARGV[1]";
while(<IN>){chomp;
	next if /^#/;
	@list=split;
	print OUT "$list[0]\t$list[1]\t$list[2]\t$list[3]\t$list[4]\n";
}
close IN;
close OUT;
