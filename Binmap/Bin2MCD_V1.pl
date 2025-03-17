######### EXAMPLE ########
##perl Bin2MCD.pl -n 12  rils_file.list  rils_traits
########################## !!!!!   For AIO-seq Pipeline  !!!20200715  #####################
use Getopt::Std;
use vars qw($opt_n $opt_o);
getopts ('n:o:');
# -n :chr number (Rice:12;  maize:10)

#####################################################################
open OUT, ">$ARGV[0].edge";
for (1..$opt_n){
	$chr_n=$_;
	open INPUT1,"<$ARGV[0]" or die "$!";
	while (<INPUT1>) {
		$filename=$_; 
		chomp($filename);
		open INPUT,"<$filename" or die "$!";
		while (<INPUT>) {
			$line=$_; 
			chomp($line);
			($chromosome,$bstart,$bend,$origin,$name)=split(/\s+/,$line);
			$bstart100=int($bstart/100000+0.5);
			$bend100=int($bend/100000+0.5);
			$number=$chromosome;
			$number=~s/chromosome//;
			$number=~s/^0//;
			if ($number == $chr_n){
				print OUT $filename."\t".$chromosome."\t".$bend."\t".$origin."\t".$bend100."\n";
			}
		}
		close INPUT;
	}   
}
close INPUT1;
close OUT;

#####################################################################
open IN, "<$ARGV[0].edge" or die "$!";
open OUT2, ">$ARGV[0].edge.sort";
$column=3;
$type="digital";
$n=0;
while (<IN>){
	$line = $_;
	@hang = split(/\s+/,$line);
	$lane_n=$#hang;
	for (0..($lane_n)){
		$list[$n][$_]=$hang[$_];
	}
	$n++;
}
if ($type eq "digital"){
	@list = sort{$a->[$column-1] <=> $b->[$column-1]} @list;
}else{
	@list = sort{$a->[$column-1] cmp $b->[$column-1]} @list;
}

$column=2;
$type="char";

if ($type eq "digital"){
	@list = sort{$a->[$column-1] <=> $b->[$column-1]} @list;
}else{
	@list = sort{$a->[$column-1] cmp $b->[$column-1]} @list;
}

for (0..($n-1)){
	$h=$_;
	for (0..($lane_n)){
		if ($_==0){
			print OUT2 $list[$h][$_];
		}else{
			print OUT2 "\t".$list[$h][$_];
		}
	}
	print OUT2 "\n";
}

open OUT3, ">$ARGV[0].edge.comb";

$edge_temp=$list[0][4];
print OUT3 $list[0][1]."\t".$list[0][4]."\n";
	for (1..($n-1)){
		$h=$_;
		if ($list[$h][4] ne $edge_temp){
			$edge_temp=$list[$h][4];
			print OUT3 $list[$h][1]."\t".$list[$h][4]."\n";
		}
	}

close OUT2;
close OUT3;

#####################################################################
$n=0;
open INPUT2,"<$ARGV[0].edge.comb" or die "$!";
while (<INPUT2>) {
	$line=$_; 
	chomp($line);
	$line = $_;
	@hang = split(/\s+/,$line);
	if ($hang[1]==0){
		next;
	}
	$lane_n=$#hang;
	for (0..($lane_n)){
		$edge[$n][$_]=$hang[$_];
	}
	$n++;
}
close INPUT2;

$chrom_temp="";
open OUT4, ">$ARGV[0].bin.mark";
for (0..($n-1)){
	$h=$_;
	if ($edge[$h][0] ne $edge_temp){
		$edge_temp=$edge[$h][0];
		$bin_mark=$edge[$h][1]/2;
		print OUT4 $edge[$h][0]."\t".$bin_mark."\n";
		next;
	}
	$bin_mark=($edge[$h][1]+$edge[$h-1][1])/2;
	print OUT4 $edge[$h][0]."\t".$bin_mark."\n";
}
close OUT4;

#####################################################################
$n=0;
@edge=();
@ril=();
open OUT5, ">$ARGV[0].map";
open INPUT3,"<$ARGV[0].bin.mark" or die "$!";
while (<INPUT3>) {
	$line=$_; 
	chomp($line);
	$line = $_;
	@hang = split(/\s+/,$line);
	$lane_n=$#hang;
	for (0..($lane_n)){
		$edge[$n][$_]=$hang[$_];
	}
	$n++;
}
close INPUT3;

$ril[0]="chr";
$ril[1]="position";
$m=2;
open INPUT,"<$ARGV[0]" or die "$!";
while (<INPUT>) {
	$filename=$_; 
	chomp($filename);
	$sn=$filename;
	$sn=~s/\.bin$//;
	print $filename."\n";
	$ril[$m]=$sn;
	open INPUT2,"<$filename" or die "$!";
	while (<INPUT2>) {
		$line=$_; 
		chomp($line);
		($chromosome,$bstart,$bend,$origin,$name)=split(/\s+/,$line);
		$bstart100=int($bstart/100000+0.5);
		$bend100=int($bend/100000+0.5);
		if ($origin eq "parent1"){
			$mark="A";
		}elsif($origin eq "parent2"){
			$mark="B";
		}elsif($origin eq "heterozygo"){
			$mark="H";
		}
		for (0..($n-1)){
			if ($chromosome eq $edge[$_][0]){
				if ($edge[$_][1]>$bstart100 && $edge[$_][1]<=$bend100){
					$edge[$_][$m]=$mark;
				}
			}
		}
	}
	close INPUT2;
	$m++;
}
$m2=$m-2;
print "RILs: $m2\n";
close INPUT;
for (0..($m-1)){
	print OUT5 $ril[$_];
	if ($_ == $m-1){
		print OUT5 "\n";
	}else{
		print OUT5 "\t";
	}
}
for (0..($n-1)){
	$h=$_;
	for (0..($m-1)){
		if ($_==0){
			print OUT5 $edge[$h][$_];
		}else{
			print OUT5 "\t".$edge[$h][$_];
		}
	}
	print OUT5 "\n";
}

close OUT5;

################################################################
$n=0;
@edge=();
$m=1;
$chrom[0]=0;
open OUT6, ">$ARGV[0].map.mcd";
open INPUT1,"<$ARGV[0].map" or die "$!";
while (<INPUT1>) {
	$line=$_; 
	chomp($line);
	$line = $_;
	@hang = split(/\s+/,$line);
	if ($hang[1] eq "position"){
		next;
	}
	$lane_n=$#hang;
	for (0..($lane_n)){
		$edge[$n][$_]=$hang[$_];
	}
	$n++;
}
close INPUT1;
$all_chromosome=$edge[$n-1][0];
$all_chromosome=~s/chromosome//;     
$all_chromosome=~s/^0//;
$sample=($lane_n)-1; #### edited by zhaosheng in 20190619  (the primary script is $sample=($lane_n)-2;)
$chromosome="chromosome01";
$max=0;
for (0..$n-1) {
	if ($edge[$_][0] eq "$chromosome"){
		$chrom[$m]++;
	}else{
		$chromosome=$edge[$_][0];
		$chrom[$m+1]++;
		$m++;
	}
}

for (1..$all_chromosome){
	if ($chrom[$_] > $max){
		$max=$chrom[$_];
	}
}

for (1..$all_chromosome){
	$chrom[$_]=$chrom[$_]+$chrom[$_-1];
}

for (1..$all_chromosome){
	$edge[$chrom[$_-1]][1]=0;
}

for (1..$all_chromosome){
	$c=$_;
	for ($chrom[$c-1]+1..($chrom[$c]-1)) {
		$h=$_;
		$count=0;
		for (2..$lane_n) {
			if (($edge[$h][$_] ne $edge[$h-1][$_]) && ($edge[$h][$_] ne "H") && ($edge[$h-1][$_] ne "H")){
				$count=$count+2;
			}elsif (($edge[$h][$_] ne $edge[$h-1][$_]) && (($edge[$h][$_] ne "H") or ($edge[$h-1][$_] ne "H"))){
				$count=$count+1;
			}
		}
		$edge[$h][1]=$edge[$h-1][1]+(100*$count)/(2*$sample);
	}
}

for (0..$n-1){
	$edge[$_][1]=int(10*$edge[$_][1]+0.5)/10;
}

print OUT6 <<_FLAG_;
#FileID	20190626
#bychromosome
 -type	position
 -function	2
 -Units	cM
 -chromosomes	$all_chromosome
 -maximum	$max
 -named	yes
 -start
_FLAG_

$sn=1;
$number_temp=0;
for (0..$n-1){
	$number=$edge[$_][0];
	$number=~s/chromosome//;     
	$number=~s/^0//;
	if ($number ne $number_temp){
		print OUT6 " -Chromosome\tC".$number."\n";
		$number_temp=$number;
	}
	if ($sn<10){
		print OUT6 "bin_000".$sn."\t".$edge[$_][1]."\n";
	}elsif ($sn>=10 && $sn<100){
		print OUT6 "bin_00".$sn."\t".$edge[$_][1]."\n";
	}elsif ($sn>=100 && $sn<1000){
		print OUT6 "bin_0".$sn."\t".$edge[$_][1]."\n";
	}elsif ($sn>=1000){
		print OUT6 "bin_".$sn."\t".$edge[$_][1]."\n";
	}
	$sn++;
}

print OUT6 <<_FLAG2_;
-stop

---------------------------------------------------
#bycross
 -SampleSize	$sample
 -Cross	Ri1
 -traits	4
 -missingtrait	.
 -case	yes
 -TranslationTable
 AA    2     2
 Aa    1     1
 aa    0     0
 A-    12    -12
 a-    10    -10
 --    -1    -1
 -start markers
_FLAG2_

$sn=1;
for (0..$n-1){
	$c=$_;
	if ($sn<10){
		print OUT6 "bin_000".$sn;
	}elsif ($sn>=10 && $sn<100){
		print OUT6 "bin_00".$sn;
	}elsif ($sn>=100 && $sn<1000){
		print OUT6 "bin_0".$sn;
	}elsif ($sn>=1000){
		print OUT6 "bin_".$sn;
	}
	
	for (2..$#hang){
		if ($edge[$c][$_] eq "A"){
			print OUT6 "\t"."0";
		}elsif ($edge[$c][$_] eq "H"){
			print OUT6 "\t"."1";
		}elsif ($edge[$c][$_] eq "B"){
			print OUT6 "\t"."2";
		}
		if ($_ eq $#hang){
			print OUT6 "\n";
		}
	}
	$sn++;
}

print OUT6 <<_FLAG3_;
 -stop	markers
 -start	traits
_FLAG3_

$n=0;
open INPUT2,"<$ARGV[1]" or die "$!";
while (<INPUT2>) {
	$line=$_; 
	chomp($line);
	$line = $_;
	@hang = split(/\s+/,$line);
	$lane_n=$#hang;
	for (0..($lane_n)){
		$trait[$n][$_]=$hang[$_];
	}
	$n++;
}
close INPUT2;
print "Traits: ".$lane_n."\n";

for (1..$lane_n){
	$t=$_;
	print OUT6 $trait[0][$t];
	for (1..$n-1){
		print OUT6 "\t".$trait[$_][$t];
	}
	print OUT6 "\n";
}
print OUT6 <<_FLAG4_;
 -stop	traits
 -quit
 -end
_FLAG4_

$bin_n=$sn-1;
print "Bin: ".$bin_n."\n";
close OUT6;

exit; 
