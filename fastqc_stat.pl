opendir (DIR, "./") or die "can't open the directory!";
@dir = readdir DIR;
foreach $file ( sort  @dir) 
{

# 跳过不需要的文件/文件夹，留下需要的文件夹
next unless -d $file;
next if $file eq '.';
next if $file eq '..';

# 提取total reads
$total_reads=  `grep '^Total' ./$file/fastqc_data.txt`;
$total_reads=(split(/\s+/,$total_reads))[2];
# 提取%GC
$GC= `grep '%GC' ./$file/fastqc_data.txt`;
$GC=(split(/\s+/,$GC))[1];
chomp $GC;
# 提取Q20，Q30
## 读入Per sequence quality scores部分的信息，保存成哈希
open FH , "<./$file/fastqc_data.txt";
while (<FH>)
    {
    next unless /#Quality/;
    while (<FH>)
        {
        @F=split;
        $hash{$F[0]}=$F[1];
        last if />>END_MODULE/;
        }
    }
## 统计Q20，Q30
$all=0;$Q20=0;$Q30=0;
$all+=$hash{$_} foreach keys %hash;
$Q20+=$hash{$_} foreach 0..20;
$Q30+=$hash{$_} foreach 0..30;
$Q20=1-$Q20/$all;
$Q30=1-$Q30/$all;

print "$file\t$total_reads\t$GC\t$Q20\t$Q30\n";
}