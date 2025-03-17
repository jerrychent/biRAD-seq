#!/usr/bin/perl -w
#usage: perl $0 ZM73  TL-Hd1_501  TL-AD2_701
my (@set,$F_R1,$F_R2,$R_R1,$R_R2,@FF,@FR,@RR,@RF,%hash1,%hash2,%hash3,%hash4);
	
$F_R1=$ARGV[1].".R1.fq.head.id"; #TL-Hd1_502.R1.fq.head.id
$F_R2=$ARGV[1].".R2.fq.head.id"; #TL-Hd1_502.R2.fq.head.id
$R_R1=$ARGV[2].".R1.fq.head.id"; #TL-AD2_702.R1.fq.head.id
$R_R2=$ARGV[2].".R2.fq.head.id"; #TL-AD2_702.R2.fq.head.id

$out=$ARGV[0]."-read.head";
open OUT,">./$out";
#===========================	
open IN1,$F_R1;
while(<IN1>){chomp;
	@set=split /\@/,$_;
	push(@FF,$set[1]);
}
close IN1;
open IN2,$F_R2;
while(<IN2>){chomp;
	@set=split /\@/,$_;
	push (@FF,$set[1]);
}
close IN2;
for(@FF){
	$hash1{$_}++;
}
for (keys %hash1){
	if($hash1{$_}==2){
		print OUT "$_\n";
	}else{
		next;
	}
}
@FF=();
#========================
open IN1,$F_R1;
while(<IN1>){chomp;
	@set=split /\@/,$_;
	push (@FR,$set[1]);
}
close IN1;
open IN4,$R_R2;
while(<IN4>){chomp;
	@set=split /\@/,$_;
	push (@FR,$set[1]);
}
close IN4;
for(@FR){
	$hash2{$_}++;
}
for (keys %hash2){
	if($hash2{$_}==2){
		print OUT "$_\n";
	}else{
		next;
	}
}
@FR=();
#=====================
open IN3,$R_R1;
while(<IN3>){chomp;
	@set=split /\@/,$_;
	push (@RR,$set[1]);
}
close IN3;
open IN4,$R_R2;
while(<IN4>){chomp;
	@set=split /\@/,$_;
	push (@RR,$set[1]);
}
close IN4;
for(@RR){
	$hash3{$_}++;
}
for (keys %hash3){
	if($hash3{$_}==2){
		print OUT "$_\n";
	}else{
		next;
	}
}
@RR=();
#=====================
open IN3,$R_R1;
while(<IN3>){chomp;
	 @set=split /\@/,$_;
	push (@RF,$set[1]);
}
close IN3;
open IN2,$F_R2;
while(<IN2>){chomp;
	@set=split /\@/,$_;
	push (@RF,$set[1]);
}
close IN2;
for(@RF){
	$hash4{$_}++;
}
for (keys %hash4){
	if($hash4{$_}==2){
		print OUT "$_\n";
	}else{
		next;
	}
}
@RF=();
#====================
close OUT;	
