#!/usr/bin/perl -w
#usage: perl $0 Offspring_SEG-map_BinFile_filter.list.map  Offspring_SEG-map_BinFile_filter.list.map.ggplot2.input  BinMap_Rscript.R
#从Offspring_SEG-map_BinFile_filter.list.map文件获得ggplot2输入文件，用于制作binmap图;同时生成R语言作图所需的R脚本BinMap_Rscript.R。
#注意！要求Offspring_SEG-map_BinFile_filter.list.map文件内仅含三种基因型（如亲本基因型A/0、杂合基因型H/1,以及第二个亲本基因型B/2），不允许第四种基因型存在！
#
#
#
my($space,@list,$len,%hash,$chr_num,$sum,$chr,$newStart,$newEnd,%start,%end);
$space=1.5;###染色体间距，此处设置为1.5Mb。可根据最终出图的间距，随时再来调整此值！
open IN,$ARGV[0];
while(<IN>){chomp;
	next if $.==1;
	@list=split /\s+/,$_,3;
	$len=$list[1]/10;     ##将100Kb坐标转换为Mb坐标
	$hash{$list[0]}=$len;
}
close IN;
$chr_num=keys %hash;
$sum=0;
####以下脚本是获得每条染色体新的起止坐标（因染色体间隙存在，后面染色体起始坐标需累加前面染色体长度及间隙长度！）======================
for(1..$chr_num){
	$chr='chromosome'.sprintf("%02d",$_);
	$newStart=($_-1)*$space+$sum;
	$newEnd=$newStart+$hash{$chr};
	$start{$chr}=$newStart;
	$end{$chr}=$newEnd;
	$sum+=$hash{$chr};
	print "$chr\t$start{$chr}\t$end{$chr}\n";
}
#####=================================================================================================================================
my($chrom,$startpos,@head,$id,$begin,$stop);
open IN,$ARGV[0];
open OUT1,">./$ARGV[1]";
print OUT1 "ID\tsample\tChr\tGenotype\tNewStart\tNewEnd\n";
$chrom="chromosome01";
$startpos=0;
while(<IN>){chomp;
	if($.==1){
		@head=split;
	}else{
		@row=split;
		if($chrom eq $row[0]){
			for(2..$#row){
				$id=$_-1;
				$begin=$start{$row[0]}+$startpos*0.1;
				$stop=$start{$row[0]}+$row[1]*0.1;
				print OUT1 "$id\t$head[$_]\t$row[0]\t$row[$_]\t$begin\t$stop\n";
			}
			$startpos=$row[1];
		}else{
			$startpos=0;
			$chrom=$row[0];
			for(2..$#row){
				$id=$_-1;
				$begin=$start{$row[0]}+$startpos*0.1;
				$stop=$start{$row[0]}+$row[1]*0.1;
				print OUT1 "$id\t$head[$_]\t$row[0]\t$row[$_]\t$begin\t$stop\n";
			}
			$startpos=$row[1];
		}
	}
}
close IN;
close OUT1;
###===============================以下脚本用于产生ggplot2作图脚本====================================================
my($color,$sam,@pop,$population,$y_centerPos,@y_coord,$y_axis);
$color='"#FF0000","#0000FF","#CCCCCC"'; ## ABH三种基因型的颜色。可 DIY ~
open OUT2,">./$ARGV[2]";
for(2..$#head){
	$sam='"'.$head[$_].'"';
	push(@pop,$sam);
}
$population=join(",",@pop);
for(1..$#head-1){
	$y_centerPos=$_+0.4; #bin高度是1-1.8区间，则样品名称标注在y=1.4位置处。
	push(@y_coord,$y_centerPos);
}
$y_axis=join(",",@y_coord);
print OUT2 "library(ggplot2)\ndata0=read.table(\"C:/Users/lenovo/Desktop/Offspring_SEG-map_BinFile_filter.list.map.ggplot2.input\",header=TRUE)\ntail(data0)\n";
print OUT2 "ggplot(data0,aes(xmin=NewStart, xmax=NewEnd, ymin=ID, ymax=ID+0.8,fill=Genotype)) +geom_rect()+theme_classic()+scale_x_continuous(expand=c(0.01,0.01),breaks=NULL)+scale_y_continuous(expand=c(0.0,0.0), breaks=c($y_axis), labels=c($population))+scale_fill_manual(values=c($color))+";
for(1..$chr_num){
	$chr='chromosome'.sprintf("%02d",$_);
	if($_%2==1){
		print OUT2 "annotate(\"rect\",xmin=$start{$chr},xmax=$end{$chr},ymin=0,ymax=0.5,color=\"#4C4C4C\",fill=\"#4C4C4C\")+\n";
	}else{
		print OUT2 "annotate(\"rect\",xmin=$start{$chr},xmax=$end{$chr},ymin=0,ymax=0.5,color=\"#999999\",fill=\"#999999\")+\n";
	}
}
print OUT2 "theme(axis.line.x = element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank())+\nlabs(x='')\n";
print OUT2 "#+annotate(\"segment\", x = 10, xend = 10, y = 0, yend = 73.8,colour = \"black\", size=1, alpha=1) ##添加竖线段";
close OUT2;
