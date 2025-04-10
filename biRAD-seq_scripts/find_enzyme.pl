#!/usr/bin/perl -w
use strict;

die "UASE:\tperl $0 cmap min\n" unless @ARGV==3;

open (IN,"$ARGV[0]")or die $!;
my %info;
my ($chr,$weizhi);
my $length;
my $total_segments = 0;      # 总片段数
my $threshold_segments = 0;  # 阈值条件下的片段数
my $total_length = 0;       # 总片段长度
my $threshold_length = 0;   # 阈值条件下的片段长度

while(<IN>){
	chomp;
	next if /^#/;
	($chr,$weizhi)=(split/\t/,$_)[0,5];
	$info{$chr}.="$weizhi\t";
	$total_segments++;  # 每读入一个片段，总片段数加一
	
}

my $num;
my @xinxi;

open (OUT,">$ARGV[0]_enzyme_$ARGV[1]-$ARGV[2]_info");
my @key=sort {$a<=>$b} keys %info;
foreach $chr(@key){
	@xinxi=(split/\t/,$info{$chr});
	for(my $i=0;$i<@xinxi-1;$i++){
		$length=$xinxi[$i+1]-$xinxi[$i];
		if ($length>=$ARGV[1] && $length<=$ARGV[2]){   # 阈值片段长度介于命令行第二个参数ARGV[1]和第三个参数ARGV[2]之间
			$num++;
			$threshold_segments++;  # 如果片段长度满足阈值条件，阈值片段数加一
			$threshold_length += $length;  # 累加阈值片段长度
			print OUT "$chr\t$xinxi[$i]\t$xinxi[$i+1]\t$length\n";
		} 
		$total_length += $length;  # 计算总片段长度
	}
}	


# 计算比例并输出
my $ratio_segments = ($threshold_segments / $total_segments) * 100;
my $ratio_length = ($threshold_length / $total_length) * 100;

# printf "Sample: $ARGV[0]\n";
# printf "Threshold Segments: %d\n", $threshold_segments;
# printf "Total Segments: %d\n", $total_segments;
# printf "Segments Ratio: %.2f%%\n", $ratio_segments;
# printf "Threshold Length: %d\n", $threshold_length;
# printf "Total Length: %d\n", $total_length;
# printf "Length Ratio: %.2f%%\n", $ratio_length;

# 输出结果，仅包含Sample和Length Ratio，以tab分隔
printf "$ARGV[0]\t%.2f%%\n", $ratio_length;