#!/usr/bin/perl -w
use strict;
#usage: perl $0 -i /path/to/Offspring_SEG-map_BinFile_CrossoverCount.out -m 200 -o /path/to/output.file
#purpose: 根据设置的允许的单株最大重组断点数目进行过滤，保留低于此值得单株。
use Getopt::Std;
use vars qw($opt_i $opt_m $opt_o);
getopts ('i:m:o:');
my(@list);
open IN,$opt_i;
open OUT,">$opt_o";
while(<IN>){chomp;
	next if /^#/;
	@list=split;
	next if ($list[$#list])>$opt_m;
	print OUT "$list[0]\n";
}
close IN;
close OUT;
