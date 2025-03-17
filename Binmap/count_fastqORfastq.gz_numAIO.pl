#!/usr/bin/perl -w
use strict;
#perl  $0  input.list  output
die "this script is used to count the number of fastq.gz in a file_list\nUSAGE:\n\tperl $0 list\n" if @ARGV==0;

my @sample;
my %num;
my @F;
my $sample;
open (LIST,"$ARGV[0]");
while (<LIST>){
	chomp;
	push (@sample,$_);
}
foreach $sample (@sample){
	if ($sample=~/gz$/){
		open (SAMPLE,"gzip -dc $sample|") or die $!;
		}
	else{
		open (SAMPLE,"$sample")or die $!;
		}
	$num{$sample}=0;
	while (<SAMPLE>){
		chomp;
		if ($.%4==2){
		s/\r//g;
		$num{$sample}+=length$_;
	}
    }
}
open(OUT,">$ARGV[1]");
foreach $sample (@sample){
	print OUT "$sample\t$num{$sample}\n";
}
			
