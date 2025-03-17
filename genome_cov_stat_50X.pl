#!/usr/bin/perl -w
#usage: perl $0 qualimap_result.list >qualimap_result.list.stat
use strict;
my(@list,$std_cov,$name,$reads_maped,$dup,$insert,$meanMQ,$mean_cov,@cov_gt);
open IN1,"<$ARGV[0]";
while(<IN1>){chomp;
    push(@list,$_);
}
close IN1;
print "Sample_ID\treads_maped(%)\tdup_rate(%)\tmeanMQ\tmean_covData(X)\tstd_covData(X)";
for my $i (1..50) {
    push(@cov_gt, "cov_gt${i}");
    print "\t${cov_gt[-1]}(%)";
}
print "\n";
my $index = 1;
foreach(@list){
    $name = $index;
    $index++;
    open IN2,"<$_";
    if(/\d \d/){print "$_\n";}
    while(<IN2>){chomp;
        if(/number of mapped reads.*(\d+,\d+,\d+) \((\d+(\.\d+)?)%\)/){$reads_maped=$2;}
        if(/duplication rate .* (\d+\.\d+)%/){$dup=$1;}
        if(/mean mapping quality .* (\d+\.\d+)/){$meanMQ=$1;}
        if(/mean coverageData.* (\d+\.\d+)X/){$mean_cov=$1;}
        if(/std coverageData.* (\d+\.\d+)X/){$std_cov=$1;}
        for my $i (1..50) {
            if(/.* (\d+(?:\.\d+)?)%.*>= ${i}X/){$cov_gt[$i-1]=$1;}
        }
    }
    print "${name}\t${reads_maped}\t$dup\t$meanMQ\t${mean_cov}\t${std_cov}";
    for my $i (0..49) {
        print "\t${cov_gt[$i]}";
    }
    print "\n";
    close IN2;
}
