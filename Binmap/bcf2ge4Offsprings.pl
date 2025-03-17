#!/usr/bin/perl -w
use Getopt::Std;
use vars qw($opt_Q $opt_q);
getopts ('Q:q:');
my(@row,%hash,$chr_pos,@line,@list,@sample,@p2,@sample_GT,@p2_GT,$dp,$label,$dp_format,$gt_format,$sample_snp,$chro,$chromo,$chromosome,$genotype);
while(<>){chomp;
	next if(/^##/);
	@line=split;
	my $chr_positon=$line[0].'_'.$line[1]; ## Gm01_394  Chr1_4654566
	if(/^#/){
		print "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$line[9]\t$line[5]\n";
		next;
	}else{
		$chro=$2 if ($line[0]=~/(.*?)(\d+)/);
		$chromo=sprintf("%02d",$chro);
		$chromosome="chromosome".$chromo;
		next if($line[5]< $opt_Q); ###根据第6列对应的QUAL值进行过滤:QUAL lower than 40 will be discarded.
		my $mq=$1 if($line[7]=~/;MQ=(\d+)/);
		next if ($mq< $opt_q);      #### MQ in INFO lower than 40 will be discard.
		my @format=split /:/,$line[8]; #拆分第9列FORMAT:   GT:DP:SP:DP4
		for $label(0..$#format){
			if($format[$label] eq 'GT'){
				$gt_format=$label;
			}else{
				next;
			}
		}
		@sample=split /:/,$line[9]; #第10列对应的第一个亲本  #1/1:41,3,0:1:0:0,0,1,0
		if($sample[$gt_format] eq '0/0'){
			$sample_snp=$line[3];
		}elsif($sample[$gt_format] eq '1/1'){	
			$sample_snp=$line[4];
		}else{
			next;
			#$sample_snp=$line[3].$line[4];
		}
		print "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$sample_snp\t$line[5]\n";
	}
}

