use strict;
use warnings;
use File::Basename;

# 合并多个样本*.exon_level_cnv_result.xls结果并输出到一个文件中

my ($list, $outfile) = @ARGV;
my @files;
open IN, "$list" or die;
while (<IN>){
    chomp;
    next if /^$/;
    push @files, $_;
}
close IN;

open O, ">$outfile" or die;
print O "sample\tgene\tchr\tstart\tend\tcopynumber\tcnv_state\ttranscript\texon\n";

for my $f (@files){
    my $name = (split /\./, basename $f)[0];
    my $line = (split /\s/, `wc -l $f`)[0];
    if ($line == 1){
        print O "$name\tNegative\t-\t-\t-\t-\t-\t-\t-\n";
    }else{
        open IN, "$f" or die;
        <IN>;
        while (<IN>){
            chomp;
            print O "$_\n";
        }
        close IN;
    }
}
close O;