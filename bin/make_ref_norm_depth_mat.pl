use strict;
use warnings;
use File::Basename;

# ref样本深度矩阵,每一行为一个target
# 输入文件:每一行一个*.normalized_depth.txt

my ($norm_depth_list, $outfile) = @ARGV;

my @file;
my @sample;
my %depth_mat;

open IN, "$norm_depth_list" or die;
while (<IN>){
    chomp;
    next if /^$/;
    my $name = (split /\./, basename $_)[0];
    push @file, $_;
    push @sample, $name;

    my $file = $_;
    open DEPTH, "$file" or die;
    <DEPTH>;
    while (<DEPTH>){
        chomp;
        my @arr = split /\t/;
        my $target = "$arr[0]\:$arr[1]\-$arr[2]"; # chr:start-end
        push @{$depth_mat{$target}}, $arr[-1];
    }
    close DEPTH;
}
close IN;

# target
my $file1 = $file[0];
my @target;
open IN, "$file1" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $target = "$arr[0]\:$arr[1]\-$arr[2]\t$arr[3]"; # chr:start-end/gene
    push @target, $target;
}
close IN;


# 表头
open O, ">$outfile" or die;
print O "target\tgene";
for (@sample){
    print O "\t$_";
}
print O "\n";

for my $t (@target){
    my @val = split /\t/, $t;
    print O "$val[0]\t$val[1]";
    my $depth_aref = $depth_mat{$val[0]};
    for my $v (@{$depth_aref}){
        print O "\t$v";
    }
    print O "\n";
}
close O;

