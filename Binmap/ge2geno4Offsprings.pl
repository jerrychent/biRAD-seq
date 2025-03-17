#!/usr/bin/perl -w
##use strict;
use Getopt::Std;
use vars qw($opt_P $opt_I);
getopts ('P:I:');
my(@row,%hash,$chr_pos,@line,$chro,$chromo,$chromosome,$genotype);
open IN0,$opt_P;
while(<IN0>){chomp;
	next if /^#/;
	@row=split;
	$chr_pos=$row[0].'_'.$row[1]; ## Gm01_394
	$hash{$chr_pos}=0;
	my $p1=$chr_pos.'P1';
	my $p2=$chr_pos.'P2';
	$hash{$p1}=$row[4]; #第一个亲本基因型
	$hash{$p2}=$row[5]; #第二个亲本基因型
	}
close IN0;
open IN1,$opt_I;
while(<IN1>){chomp;
	@line=split;
	my $chr_positon=$line[0].'_'.$line[1]; ## Gm01_394  Chr1_4654566
	my $p3=$chr_positon.'P1';
	my $p4=$chr_positon.'P2';
	if(/^#/){
		print "#P1P2H\tChr_position\tChromosome\tPOS\tPOS\t$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\tP1\tP2\t$line[5]\n";
	}elsif(defined $hash{$chr_positon}){
		$chro=$2 if ($line[0]=~/(.*?)(\d+)/); ## 1
		$chromo=sprintf("%02d",$chro);    ## 01
		$chromosome="chromosome".$chromo;
		if($line[4] eq $hash{$p3}){
			$genotype="P1";
		}elsif($line[4] eq $hash{$p4}){	
			$genotype="P2";
		}else{
			next;
		}
		print "$genotype\t$chr_positon\t$chromosome\t$line[1]\t$line[1]\t$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$hash{$p3}\t$hash{$p4}\t$line[5]\n";
	}else{
		next;
	}
}
close IN1;
