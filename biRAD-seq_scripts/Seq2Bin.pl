######### EXAMPLE ########
##perl Seq2Bin.pl Lane1.AAT.PE.fastq.rlt jap_v4_length_list
##########################
use strict;
#use warnings;
my @hang;
my @array;
my @chrom;
my @bin;
my @len_detail;
my @temp;

my $parent1=0;
my $parent2=0;
my $win_size=0;
my $head=0;
my $head_edge=0;
my $max_n=0;
my $hetero_12="";
my $hetero_key="";
my $hetero_start="";
my $hetero_end="";
my $line="";
my $line1="";
my $chromosome1="";
my $cstart="";
my $cend="";
my $number1="";
my $chrom_len="";
my $lane_n=0;
my $bstart="";
my $bend="";
my $origin="";
my $number="";
my $c=0;
my $h=0;
my $len=0;
my $round=0;
my $name="";
#####################################################
#Judge the edges of every bin.#
#####################################################
open INPUT,"<$ARGV[0]" or die "$!";  ###  P2     HWUSI-EAS731_71:7:33:76:954#0/2      chromosome02        28560263         28560263

my $n=0;
my $m=1;
$chrom[0]=0;
my $chromosome="chromosome01";
while (<INPUT>) {
	$line = $_;
	chomp($line);
	@hang = split(/\s+/,$line);
	$lane_n=$#hang;
	for (0..($lane_n)){				###遍历此行的每一个元素的索引。
		$array[$n][$_]=$hang[$_];  ###第0行（$n）的第0个元素（$_）对应的值是 $hang[$_] .
	}
	$n++;  ### 下一行（即下一个SNP的信息）
}
if ($n<5000){						###如果此单株SNP个数小于5000个
	$win_size=15;
}elsif($n>=5000 && $n<10000){
	$win_size=25;
}elsif($n>=10000 && $n<20000){
	$win_size=49;
}elsif($n>=20000 && $n<50000){
	$win_size=99;
}elsif($n>=50000){
	$win_size=199;
}
my $half_win=int($win_size/2);  ### 7
my $win_size_dom=int($win_size*2/3);  #### 15*2/3=10    75%基因型判定标准。
open INPUT2,"<$ARGV[1]" or die "$!";
while (<INPUT2>) {                     ###遍历所有染色体
	$line1=$_; 
	chomp($line1);
	($chromosome1,$chrom_len)=split(/\s+/,$line1);
	$number1=$chromosome1;
	$number1=~s/chromosome//;
	$number1=~s/^0//;    			###染色体编号变成 1 2 3 4 5 6 7 8 9 10 11 。。。
	$len_detail[$number1]=$chrom_len;
}
my $all_chromosome=$number1;        ###最后的染色体编号是number1
close INPUT2;

for (0..$n) {                         ###遍历所有行（即判断此单株所有SNP信息）
	if ($array[$_][2] eq "$chromosome"){   ###line48，即1号染色体
		$chrom[$m]++;                      ###统计每条染色体上的SNP总个数
	}else{
		$chromosome=$array[$_][2];  ###染色体号变为下一条染色体
		$chrom[$m+1]++;
		$m++;
	}
}
for (1..$all_chromosome){
	$chrom[$_]=$chrom[$_]+$chrom[$_-1];   #### 
}

for (1..$all_chromosome){      ###遍历所有染色体
	$c=$_;
	for (($chrom[$c-1]+$half_win)..($chrom[$c]-$half_win-1)) {  ## 7..（100-7-1）
		my $win_start=$_;										## 7
		for ((($win_start)-$half_win)..($win_start+$half_win)){  ### 0..14   分别统计此窗口内P1及P2基因型的总数
			if ($array[$_][0] eq "P1"){
				$parent1++;
			}elsif($array[$_][0] eq "P2"){
				$parent2++;
			}
		}
		
		if ($parent1>$parent2){
			$array[$win_start][5]="parent1";   ### $array[7][5]="parent1"  此行多了一个元素 “parent1”
		}elsif ($parent1<$parent2){
			$array[$win_start][5]="parent2";
		}
		if ($parent1>$win_size_dom){			###如果此窗口内P1 SNP数目大于10个
			$array[$win_start][6]="parent1";     ## $array[7][6]="parent1"  此行后又多了一个元素 “parent1”
			$array[$win_start][7]=$parent1.":".$parent2; ##再在行尾加上P1和P2 SNP个数的比值 ？：？
		}elsif ($parent1<($win_size-$win_size_dom)){  ###如果P1 SNP个数小于15-10=5个
			$array[$win_start][6]="parent2";
			$array[$win_start][7]=$parent1.":".$parent2;
		}else {
			$array[$win_start][6]="hetero";
			$array[$win_start][7]=$parent1.":".$parent2;
		}
### Get the half-win of every chromosome's fisrt bin.
		if ($win_start eq $chrom[$c-1]+$half_win){
RRR1:			for ($chrom[$c-1]..$win_start-1){
				if (($array[$win_start][5] eq "parent1") && ($array[$_][0] eq "P1")){
					$head=$_;
					for ($head..$head+$half_win-1){
						$array[$_][5]="parent1";
					}
					last RRR1;
				}
				if (($array[$win_start][5] eq "parent2") && ($array[$_][0] eq "P2")){
					$head=$_;
					for ($head..$head+$half_win-1){
						$array[$_][5]="parent2";
					}
					last RRR1;
				}
				if (($array[$win_start][5] eq "parent1") && ($array[$_][0] eq "P2")){
					$head=$_;
					for ($head+1..$head+$half_win){
						if($array[$_][0] ne $array[$_-1][0]){
							$head_edge=$_;
							if ($head_edge==$head+1){
								for ($head..$head+6){
									$array[$_][5]="parent1";
								}
								last RRR1;
							}
							for ($head..$head_edge-1){
								$array[$_][5]="parent2";
							}
							for ($head_edge..$head+6){
								$array[$_][5]="parent1";
							}
							last RRR1;
						}
					}
				}
				if (($array[$win_start][5] eq "parent2") && ($array[$_][0] eq "P1")){
					$head=$_;
					for ($head+1..$head+$half_win){
						if($array[$_][0] ne $array[$_-1][0]){
							$head_edge=$_;
							if ($head_edge==$head+1){
								for ($head..$head+$half_win-1){
									$array[$_][5]="parent2";
								}
								last RRR1;
							}
							for ($head..$head_edge-1){
								$array[$_][5]="parent1";
							}
							for ($head_edge..$head+$half_win-1){
								$array[$_][5]="parent2";
							}
							last RRR1;
						}
					}
				}
			}
		}
		
### Get the half-win of every chromosome's last bin.
		my $i=0;
		my $j=0;		
		if ($win_start eq $chrom[$c]-$half_win-1){
RRR2:			for($i=$chrom[$c]-1;$i>=$win_start+1;$i--) {
				$head=$i;
				if (($array[$win_start][5] eq "parent1") && ($array[$head][0] eq "P1")){
					for ($head-$half_win+1..$head){
						$array[$_][5]="parent1";
					}
					last RRR2;
				}
				if (($array[$win_start][5] eq "parent2") && ($array[$head][0] eq "P2")){
					for ($head-$half_win+1..$head){
						$array[$_][5]="parent2";
					}
					last RRR2;
				}
				if (($array[$win_start][5] eq "parent1") && ($array[$head][0] eq "P2")){
					for($j=$head-1;$j>=$head-$half_win-1;$j--) {
						if($array[$j][0] ne $array[$j+1][0]){
							$head_edge=$j;
							if ($head_edge==$head-1){
								for ($head-$half_win+1..$head){
									$array[$_][5]="parent1";
								}
								last RRR2;
							}
							for ($head_edge+1..$head){
								$array[$_][5]="parent2";
							}
							for ($head-$half_win+1..$head_edge){
								$array[$_][5]="parent1";
							}
							last RRR2;
						}
					}
				}
				if (($array[$win_start][5] eq "parent2") && ($array[$head][0] eq "P1")){
					for($j=$head-1;$j>=$head-$half_win-1;$j--) {
						if($array[$j][0] ne $array[$j+1][0]){
							$head_edge=$j;
							if ($head_edge==$head-1){
								for ($head-$half_win+1..$head){
									$array[$_][5]="parent2";
								}
								last RRR2;
							}
							for ($head_edge+1..$head){
								$array[$_][5]="parent1";
							}
							for ($head-$half_win+1..$head_edge){
								$array[$_][5]="parent2";
							}
							last RRR2;
						}
					}
				}
			}
		}
		$parent1=0;
		$parent2=0;
	}
}

### Adjust the edges.
my $key1=1;
my $key2=1;
my $max=0;
for (0..$chrom[$all_chromosome]) {
	my $edge=$_;
	if (($array[$edge][5] ne "") && ($array[$edge+1][5] ne "") && ($array[$edge][5] ne $array[$edge+1][5]) && ($array[$edge][2] eq $array[$edge+1][2]) && ($array[$edge][6] ne "")){
		for(1..10){
			my $temp=$_-5;
			for(1..4){
				if ($array[$edge-$_+$temp][0] eq $array[$edge+$temp][0]){
					$key1++;
				}elsif($array[$edge-$_+$temp][0] ne $array[$edge+$temp][0]){
					last;
				}
			}
			for(2..5){
				if ($array[$edge+$_+$temp][0] eq $array[$edge+1+$temp][0]){
					$key2++;
				}elsif($array[$edge+$_+$temp][0] ne $array[$edge+1+$temp][0]){
					last;
				}
			}
			$array[$edge+$temp][8]=$key1+$key2;
			if ($array[$edge+$temp][8]>$max){
				$max=$array[$edge+$temp][8];
				$max_n=$edge+$temp;
			}
			$key1=1;
			$key2=1;
		}
		if (($array[$edge][8]>=5) && ($array[$max_n+1][6] ne "")){
			$array[$max_n][9]=$array[$edge][5];
			$array[$max_n+1][9]=$array[$edge+1][5];
			$max=0;
		}elsif(($array[$edge][8]>=5) && ($array[$max_n+1][6] eq "")){
			$array[$edge][9]=$array[$edge][5];
			$array[$edge+1][9]=$array[$edge+1][5];
			$max=0;
		}elsif($array[$edge][8]<5){
			$array[$edge][9]=$array[$edge][5];
			$array[$edge+1][9]=$array[$edge+1][5];
			$max=0;
		}
	}
	if (($array[$edge][5] ne "") && ($array[$edge+1][5] ne "") && ($array[$edge][5] ne $array[$edge+1][5]) && ($array[$edge][2] eq $array[$edge+1][2]) && ($array[$edge][6] eq "")){
		$array[$edge][9]=$array[$edge][5];
		$array[$edge+1][9]=$array[$edge+1][5];
	}
}
### Judge the heterozygous region.
RRR3: for (0..$chrom[$all_chromosome]) {
	if (($array[$_][6] eq "hetero") && ($array[$_-1][6] ne "hetero")){
		$hetero_start=$_;
		$hetero_key=0;
		$hetero_12=$array[$hetero_start][5];
		next RRR3;
	}
	if (($array[$_][6] eq "hetero") && ($array[$_-1][6] eq "hetero")){
		if ($hetero_12 eq $array[$_][5]){
			next RRR3;
		}elsif ($hetero_12 ne $array[$_][5]){
			$temp[$hetero_key]=$_;
			$hetero_key++;
			
			$hetero_12=$array[$_][5];
			next RRR3;
		}
	}
	if (($array[$_][6] ne "hetero") && ($array[$_-1][6] eq "hetero")){
		$hetero_end=$_-1;
		if ($hetero_key==2 && ($temp[1]-$temp[0]<int($half_win/2) )){
			$array[$hetero_start-1][9]=$array[$hetero_start-1][6];
			$array[$hetero_end+1][9]=$array[$hetero_end+1][6];
			$array[$hetero_start][9]="heterozygo";
			$array[$hetero_end][9]="heterozygo";
			for ($hetero_start+1..$hetero_end-1){
				$array[$_][9]="";
			}
		}
		if ($hetero_key>2){
			$array[$hetero_start-1][9]=$array[$hetero_start-1][6];
			$array[$hetero_end+1][9]=$array[$hetero_end+1][6];
			$array[$hetero_start][9]="heterozygo";
			$array[$hetero_end][9]="heterozygo";
			for ($hetero_start+1..$hetero_end-1){
				$array[$_][9]="";
			}
		}
		$hetero_key=0;
		$hetero_start="";
		$hetero_end="";
		next RRR3;
	}
}

open OUT1, ">$ARGV[0].win$win_size.edge";

for (0..$chrom[$all_chromosome]){
	$h=$_;
	for (0..9){
		if ($_==0){
			print OUT1 $array[$h][$_];
		}else{
			print OUT1 "\t".$array[$h][$_];
		}
	}
	print OUT1 "\n";
}			
			
close INPUT;
close OUT1;

#####################################################
#Get every bin based on the edges.#
#####################################################

open IN2,"<$ARGV[0].win$win_size.edge" or die "$!";
my $filename=$ARGV[0].".win$win_size.edge";
my ($prefix,$cenfix,$suffix)=split(/\./,$filename);
open OUT2, ">$prefix.$cenfix.bin";
$n=0;
$m=1;
@array=();
@chrom=();
$chrom[0]=0;
@bin=();
my $ij="";
my $sn=0;
$chromosome="chromosome01";
while (<IN2>) {
	$line = $_;
	chomp($line);
	@hang = split(/\t/,$line);
	$lane_n=9;
	for (0..($lane_n)){
		$array[$n][$_]=$hang[$_];
	}
	$n++;
}
close IN2;
for (0..$n) {
	if ($array[$_][2] eq "$chromosome"){
		$chrom[$m]++;
	}else{
		$chromosome=$array[$_][2];
		$chrom[$m+1]++;
		$m++;
	}
}
for (1..$all_chromosome){
	$chrom[$_]=$chrom[$_]+$chrom[$_-1];
}
for (1..$all_chromosome){
	$c=$_;
	
	for ($chrom[$c-1]..$chrom[$c]-1) {
		if ($_==$chrom[$c-1]){
			$bin[$sn][0]=$array[$_][2];
			$bin[$sn][1]=1;
			$ij=$array[$_][5];
		}
		if (($array[$_][9] ne "") && ($array[$_+1][9] ne "") && ($array[$_][9] ne $array[$_+1][9])){
			$bin[$sn][2]=int(($array[$_][3]+$array[$_+1][3])/2);
			$bin[$sn][3]=$array[$_][9];
			$bin[$sn][4]=$array[$_][1];
			$bin[$sn][5]=$bin[$sn][2]-$bin[$sn][1]+1;
			$sn++;
			$bin[$sn][0]=$array[$_][2];
			$bin[$sn][1]=$bin[$sn-1][2]+1;
			$ij=$array[$_+1][9];
		}
		if ($_==$chrom[$c]-1){
			$bin[$sn][2]=$len_detail[$c];
			$bin[$sn][3]=$ij;
			$bin[$sn][4]="chr_end";
			$bin[$sn][5]=$bin[$sn][2]-$bin[$sn][1]+1;
			$sn++;
		}
	}
}

for (0..($sn-1)){
	$h=$_;
	for (0..5){
		if ($_==0){
			print OUT2 $bin[$h][$_];
		}else{
			print OUT2 "\t".$bin[$h][$_];
		}
	}
	print OUT2 "\n";
}
close OUT2;

#####################################################
# Filtering bins which are smaller than 300K. #
#####################################################
$n=0;
for(1..100){
	$round=$_;
	open IN4,"<$prefix.$cenfix.bin" or die "$!";
	$n=0;
	while (<IN4>) {
		$line = $_;
		@hang = split(/\s+/,$line);
		$lane_n=4;
		for (0..($lane_n)){
			$bin[$n][$_]=$hang[$_];
		}
		$n++;
	}
	close IN4;
	open OUT3, ">$prefix.$cenfix.bin";
	my $skip=0;
	for (0..($n-1)){
		$h=$_+$skip;
		if ($h<$n-1) {
			if (($bin[$h+1][2]-$bin[$h+1][1]+1)>=300000){
				for (0..4){
					if ($_==0){
						print OUT3 $bin[$h][$_];
					}else{
						print OUT3 "\t".$bin[$h][$_];
					}
				}
				$len=$bin[$h][2]-$bin[$h][1]+1;
				print OUT3 "\t".$len;
				print OUT3 "\n";
			}elsif((($bin[$h+1][2]-$bin[$h+1][1]+1)<300000) && ($bin[$h+2][3] eq $bin[$h][3])){
				print OUT3 $bin[$h][0]."\t";
				print OUT3 $bin[$h][1]."\t";
				if (($bin[$h+1][4] ne "chr_end") && ($bin[$h+1][1] ne "1")){
					print OUT3 $bin[$h+2][2]."\t";
					print OUT3 $bin[$h+2][3]."\t";
					print OUT3 $bin[$h+2][4]."\t";
					$len=$bin[$h+2][2]-$bin[$h][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
					$skip=$skip+2;
				}elsif ($bin[$h+1][1] eq "1"){
					print OUT3 $bin[$h][2]."\t";
					print OUT3 $bin[$h][3]."\t";
					print OUT3 $bin[$h][4]."\t";
					$len=$bin[$h][2]-$bin[$h][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
				}elsif ($bin[$h+1][4] eq "chr_end"){
					print OUT3 $bin[$h][2]."\t";
					print OUT3 $bin[$h][3]."\t";
					print OUT3 $bin[$h][4]."\t";
					$len=$bin[$h][2]-$bin[$h][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
					
					print OUT3 $bin[$h+1][0]."\t";
					print OUT3 $bin[$h+1][1]."\t";
					print OUT3 $bin[$h+1][2]."\t";
					print OUT3 $bin[$h+1][3]."\t";
					print OUT3 $bin[$h+1][4]."\t";
					$len=$bin[$h+1][2]-$bin[$h+1][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
					$skip=$skip+1;
				}
			}elsif((($bin[$h+1][2]-$bin[$h+1][1]+1)<300000) && ($bin[$h+2][3] ne $bin[$h][3])){
				for (0..4){
					if ($_==0){
						print OUT3 $bin[$h][$_];
					}else{
						print OUT3 "\t".$bin[$h][$_];
					}
				}
				$len=$bin[$h][2]-$bin[$h][1]+1;
				print OUT3 "\t".$len;
				print OUT3 "\n";
			}
		}
		elsif ($h==$n-1) {
			for (0..4){
				if ($_==0){
					print OUT3 $bin[$h][$_];
				}else{
					print OUT3 "\t".$bin[$h][$_];
				}
			}
			$len=$bin[$h][2]-$bin[$h][1]+1;
			print OUT3 "\t".$len;
			print OUT3 "\n";
		}
	}
	close OUT3;
}

#####################################################
#Draw a figure in PNG format indicating bins and SNPs of the RIL.
#####################################################
open OUT4, ">$prefix.$cenfix.combine.png";
use GD; 
my $im = new GD::Image(5000,2000); 
my $black = $im->colorAllocate(0,0,0); 
my $white = $im->colorAllocate(255,255,255); 
my $red = $im->colorAllocate(255,0,0);
my $blue = $im->colorAllocate(0,0,255);
my $green = $im->colorAllocate(0,255,0);
my $hetero = $im->colorAllocate(255,255,0);
my $background=$im->colorAllocate(200,200,200);
$im->fill(10,10,$white); 
$im->rectangle(0,20,5000,30,$black);
$im->fill(2500,25,$black);
for (1..50){
	$im->line($_*100,5,$_*100,20,$black);
}
for (1..9){
	my $scale_temp=$_*5;
	my $scale=$scale_temp." Mb";
	$im->string(gdGiantFont,$_*500-25,50,"$scale",$black);
}

open INPUT2,"<$ARGV[1]" or die "$!";
while (<INPUT2>) {
	$line1=$_; 
	chomp($line1);
	($chromosome1,$chrom_len)=split(/\s+/,$line1);
	$number1=$chromosome1;
	$number1=~s/chromosome//;
	$number1=~s/^0//;
	$cend=int(($chrom_len/10000)+0.5);
	$im->rectangle(0,$number1*120,$cend,($number1*120+100),$background);
	$im->fill($cend*0.5,($number1*120+50),$background); 
}
close INPUT2;

open IN5,"<$prefix.$cenfix.bin" or die "$!";
while (<IN5>) {
	$line=$_; 
	chomp($line);
	($chromosome,$bstart,$bend,$origin,$name)=split(/\s+/,$line);
	$bstart=int($bstart/10000);
	$bend=int($bend/10000);
	$number=$chromosome;
	$number=~s/chromosome//;
	$number=~s/^0//;
	if ($origin eq "parent1"){
		$im->rectangle($bstart,$number*120+21,$bend,($number*120+40),$blue);
		$im->fill(($bstart+$bend)*0.5,($number*120+30),$blue); 
	}elsif ($origin eq "parent2"){
		$im->rectangle($bstart,$number*120+21,$bend,($number*120+40),$red);
		$im->fill(($bstart+$bend)*0.5,($number*120+30),$red); 
	}
	elsif ($origin eq "heterozygo"){
		$im->rectangle($bstart,$number*120+21,$bend,($number*120+40),$hetero);
		$im->fill(($bstart+$bend)*0.5,($number*120+30),$hetero); 
	}
}
close IN5;

open INPUT4,"<$ARGV[0]" or die "$!";
while (<INPUT4>) {
	$line=$_; 
	chomp($line);
	($origin,$name,$chromosome,$bstart,$bend)=split(/\s+/,$line);
	$bstart=int($bstart/10000);
	$bend=int($bend/10000);
	$number=$chromosome;
	$number=~s/chromosome//;
	$number=~s/^0//;
	if ($origin eq "P1"){
		$im->rectangle($bstart,$number*120+51,$bend,($number*120+70),$blue);
		$im->fill(($bstart+$bend)*0.5,($number*120+60),$blue); 
	}elsif ($origin eq "P2"){
		$im->rectangle($bstart,$number*120+71,$bend,($number*120+90),$red);
		$im->fill(($bstart+$bend)*0.5,($number*120+90),$red); 
	}
}
close INPUT4;
	binmode OUT4; 
	print OUT4 $im->png; 
close OUT4;
exit; 
