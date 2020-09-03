use strict;
use warnings;

# 合并多个样本*.QC.xls
my ($qc_list, $outfile) = @ARGV;

my @file;
open IN, "$qc_list" or die;
while (<IN>){
    chomp;
    next if /^$/;
    push @file, $_;
}
close IN;

# 表头
my $file1 = $file[0];
open IN, "$file1" or die;
my $h = <IN>;
close IN;

open O, ">$outfile" or die;
print O "$h";

for my $f (@file){
    print "$f\n";
    open IN, "$f" or die;
    <IN>;
    my $line = <IN>;
    close IN;
    print O "$line";
}
close O;