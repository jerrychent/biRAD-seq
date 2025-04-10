#!/usr/bin/perl -w
##usage: perl $0 DTH2_F01.R1.fq.head.id  DTH2_F01.R2.fq.head.id  DTH2_F01.R1-R2.fq.head.id.common
my (%hash0);
open IN0,$ARGV[0];
while(<IN0>){chomp;
	$hash0{$_}=0;
}
close IN0;
open IN1,$ARGV[1];
open OUT,">./$ARGV[2]";
while(<IN1>){chomp;
	if(defined $hash0{$_}){
		print OUT "$_\n";
	}else{
		next;
	}
}
close IN1;
close OUT;
