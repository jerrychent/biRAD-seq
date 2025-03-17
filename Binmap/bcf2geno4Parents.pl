#!/usr/bin/perl -w
use strict;
#usage: bcftools view parental_snp_chrall.bcf |perl $0 -U 60 -L 6 -A 30 -a 4 -B 45 -b 4 -Q 40 -q 40  >parental_snp.geno
#目的： 过滤bcf文件中低质量的SNP。
#Data:2020-07-10
#Author:dashengzhao@126.com
#============= 参数解读 =================================================================================================================================================================
use Getopt::Std;
#my %options;
use vars qw($opt_U $opt_L $opt_A $opt_a $opt_B $opt_b $opt_Q $opt_q);
getopts ('U:L:A:a:B:b:Q:q:');
#EXAMPLE:
# -U 60 意思是为避免基因组repeat区的SNP，粗略要求某一个SNP位点的reads总覆盖度（INFO列DP对应的数值），最高不超过两个亲本重测序总深度的3倍（假如双亲测序深度分别为10X和15X，则此处应该设置为75）
# -L 6  意思是为避免基因组部分区段的极偏比对，粗略要求某一个SNP位点的reads总覆盖度（INFO列DP对应的数值），最低不低于两个亲本重测序总深度的0.3倍（假如双亲测序深度分别为10X，则此处应该设置为8）。一般此值设置不应该低于6！
# -A 30  要求vcf文件第10列对应的亲本，在此SNP位点的reads覆盖深度最高不高于此亲本重测序深度的3倍；10*3=30
# -a 4   要求vcf文件第10列对应的亲本，在此SNP位点的reads覆盖深度最低不低于此亲本重测序深度的0.3倍；10*0.3=3.3 （一般此最小值大于4即可）
# -B 45  要求vcf文件第11列对应的亲本，在此SNP位点的reads覆盖深度最高不高于此亲本重测序深度的3倍；15*3=45
# -b 4   要求vcf文件第11列对应的亲本，在此SNP位点的reads覆盖深度最低不低于此亲本重测序深度的0.3倍；15*0.3=4.5 （一般此最小值大于4即可）
# -Q 40  要求vcf文件第6列QUAL值高于此处设置的值（此处为40）
# -q 40  要求vcf文件第8列INFO列中的MQ值高于此处设置的值（此处为40）
#=======================================================================================================================================================================================
my (@line,@list,@p1,@p2,@p1_GT,@p2_GT,$dp,$label,$dp_format,$gt_format,$p1_snp,$p2_snp);
while(<>){chomp;
	next if(/^##/);
	@line=split;
	@list=split /,/,$line[4];
	next if(@list>1);  #只保留二等位的SNP位点
	if(/^#/){
		print "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$line[9]\t$line[10]\t$line[5]\n"; #输出7列信息：#CHROM  POS   REF   ALT   P1_ID   P2_ID   QUAL
	}else{	next if($line[5]< $opt_Q); ###根据第6列对应的QUAL值进行过滤:QUAL lower than 40 will be discarded.
		my $mq=$1 if($line[7]=~/;MQ=(\d+)/);
		next if ($mq< $opt_q);      #### MQ in INFO lower than 40 will be discard.
		$dp=$1 if($line[7]=~/DP=(\d+);/);
		next if($dp>$opt_U or $dp<$opt_L);### the depth of all raw reads at this site should higher than 6 and lower than 60.
		my @format=split /:/,$line[8]; #拆分第9列
		for $label(0..$#format){
			if($format[$label] eq 'GT'){
				$gt_format=$label;
			}elsif($format[$label] eq 'DP'){
				$dp_format=$label;
			}else{
				next;
			}
		}
		@p1=split /:/,$line[9]; #第10列对应的第一个亲本
		@p2=split /:/,$line[10];#第11列对应的第二个亲本
		next if($p1[$gt_format] eq $p2[$gt_format]);##to discrad same genotype between two parents,即双亲无多态性的SNP标记删除。
		next if($p1[$dp_format]<$opt_a or $p1[$dp_format]>$opt_A);###对于第10列的亲本，the depth of qualified reads at this site should range from $opt_a to $opt_A.
		next if($p2[$dp_format]<$opt_b or $p2[$dp_format]>$opt_B);###对于第11列的亲本，the depth of qualified reads at this site should range from $opt_b to $opt_B.
		@p1_GT=split /\//,$p1[$gt_format];
		next unless($p1_GT[0]==$p1_GT[1]);##to discard heterozygous genotype
		@p2_GT=split /\//,$p2[$gt_format];
		next unless($p2_GT[0]==$p2_GT[1]);##to discard heterozygous genotype
		if($p1[$gt_format] eq '0/0'){
			$p1_snp=$line[3];
		}elsif($p1[$gt_format] eq '1/1'){
			$p1_snp=$line[4];
		}else{
			$p1_snp="N";
		}
		if($p2[$gt_format] eq '0/0'){
			$p2_snp=$line[3];
		}elsif($p2[$gt_format] eq '1/1'){
			$p2_snp=$line[4];
		}else{
			$p2_snp="N";
		}
		print  "$line[0]\t$line[1]\t$line[3]\t$line[4]\t$p1_snp\t$p2_snp\t$line[5]\n";  #输出7列信息：#CHROM  POS   REF   ALT   P1_ID   P2_ID   QUAL
		
	}
}
###===== have a nice day, my friend ! ===========
