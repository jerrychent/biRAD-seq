#!/usr/bin/perl -w
#changed in 2021-12-21 by ZS
#usage: perl $0 
#usage: perl $0 -i R498_reads_q30.sort.RGadd.rmdu.bam.R1.plusStrand.bam.depth  -d 3 -o find_R1_plusStrand_tags.out1
#-i: input
#(disabled: -c : the coverage depth of whole raw data.)
#-d: the depth of tags site should bigger than this value.)
#-o: output
#purpose: count tags and its coordinates on plusStrand.
my(@list,%hash1,%hash2,%hash3,@row,$name,$hou1,$all_line,$max,@hang,$qian1);
use Getopt::Std;
use vars qw($opt_i  $opt_d $opt_o);
getopts ('i:d:o:');
open IN,$opt_i;
while(<IN>){chomp;
	@list=split;
	$pos=$list[0]."-".$list[1]; #chr1-12345
	$hash1{$.}=$list[2];#line2=>dep1
	$hash2{$.}=$_;      #line2=>wholeLine
	$hash3{$pos}=$_;    #chr1-12345=>wholeLine
}
close IN;
open OUT,">./$opt_o";
#print OUT "#Chr\tPos(1-based)\tdepth(x)\tR1_match_to_+/-_Strand\n";
print OUT "#Chr\tPos(1-based)\tdepth(x)\n";
my $line=1;
$all_line=keys %hash2;
while($line<$all_line){ 
	if($hash1{$line}>=$opt_d){ #depth >= $opt_d
		@list=split /\s+/,$hash2{$line}; #Chr1    12221   10
		$hou1=$list[1]-1;                #12220
		$name=$list[0]."-".$hou1;        #Chr1-12220
		$qian1=$list[1]+1;
		$id=$list[0]."-".$qian1;        #Chr1-12222
		if(defined $hash3{$name}){       #Chr1-22220=>wholeLine exist.
			@row=split /\s+/,$hash3{$name};
		}else{
			$row[2]=0; #dep=0
			$row[0]=$list[0]; #chr=Chr1     @row="Chr1 12220 0"
		}

		if(defined $hash3{$id}){    #Chr1-22222=>wholeLine  exist.
			@hang=split /\s+/,$hash3{$id};
		}else{
			$hang[2]=0; #dep=0
			$hang[0]=$list[0]; #chr=Chr1     @hang="Chr1 12222 0"
		}

		if(($list[0] eq $row[0]) and  ($list[2]>=$row[2]+$opt_d) ){
			print OUT "$hash2{$line}\n";
			$line++;
		}elsif(($list[0] eq $hang[0]) and ($list[2]>=$hang[2]+$opt_d)){
			print OUT "$hash2{$line}\n";
			$line++;
		}else{
			$line++;
		}
	}else{
		$line++;
	}
}
close OUT;
