#!/usr/bin/perl -w
#usage: perl $0 TL-seqA1_PE140_q30.sort.bam.dep.tags.reads  TL-seqA1_PE140_q30.sort.bam.dep.tags.reads.win
my (@line,$chr,$tag,$pos,$int,$right,$win,%hash,@set,$left,$reg,$total,%haChr,%haPos,@head);
open IN0,"$ARGV[0]"; #TL-seqA1_PE140_q30.sort.bam.dep.tags.reads
while(<IN0>){chomp;
	next if /^#/;
	@line=split;
	$chr=$1 if ($line[0]=~/hr(.*)$/);
	#$chr=$line[0];
	$tag=$line[0]."_".$line[1];
	$pos=$line[1]/1000;
	$int=int $pos;
	$right=$int +1;
	$win=$chr."-".$right;
	$hash{$win}++;
	push(@{$win},$tag);
}
close IN0;
open OUT0, ">./$ARGV[1].tem";
#print OUT0 "winID\tRegion(Kb)\ttagNum\ttagCoord\n";
foreach(sort keys%hash){
	@set=split /-/,$_;
	$left=$set[1]-1;
	$reg=$left."<=Win<".$set[1];
	$total=join(";",@{$_});
	print OUT0 "$_\t$reg\t$hash{$_}\t$total\n";
	$total="";
}
close OUT0;
open IN1,"$ARGV[1].tem";
open OUT1,">./$ARGV[1]";
while(<IN1>){chomp;
	next if /^#/;
	@line=split /\s+/,$_,2;
	@head=split /-/,$line[0];
	
	$haChr{$head[0]}=0;
	$haPos{$head[0]}{$head[1]}=$_;
}
close IN1;
open OUT1,">./$ARGV[1]";
print OUT1 "#winID\tRegion(Kb)\ttagNum\ttagCoord\n";
for $chr (sort {$a <=> $b} keys %haChr){
	for $pos (sort {$a <=> $b} keys %{$haPos{$chr}}){
		print OUT1 "$haPos{$chr}{$pos}\n";
	}
}
close OUT1;
unlink "$ARGV[1].tem";
