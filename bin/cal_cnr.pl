use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw/$Bin/;
use List::Util qw(sum);
use POSIX qw(ceil);
use POSIX qw(floor);

# 2020/1/30

my ($file,$ctr_num,$ctr_dir,$outdir);

GetOptions(
    "infile:s" => \$file,           # NEED
    "ctr_num:i" => \$ctr_num,       # Optional
    "ctr_dir:s" => \$ctr_dir,       # NEED
    "outdir:s" => \$outdir,         # NEED
    ) or die "unknown args\n";

# check args
if (!-e $file){
    die "can not find *.normalized_depth.txt file\n";
}

# default value
if (not defined $ctr_num){
    $ctr_num = 70;
}

# check ctr_dir
if (!-d $ctr_dir){
    die "can not find control dir\n";
}


my $name = (split /\./, basename $file)[0];
print "auto detected sample name (from norm depth file) is: $name\n";
print "cal copy number ratio for $name\n";

# get all ref sample
my @control_files = glob "$ctr_dir/*.normalized_depth.txt";
if (@control_files == 0){
    # no ref
    die "can not find *.normalized_depth.txt file under $ctr_dir\n";
}

my $n = scalar @control_files;
print "find $n refs under ctr_dir\n";

if ($n < $ctr_num){
    print "WARNING: there are only $n control files.\n";
    print "all these control files will be used as refs\n";
    $ctr_num = $n;
}

# cal abs diff
my %avgDIFF;
for my $nc (@control_files){
    my $avgdiff = &cal_cov_diff($file, $nc);
    $avgDIFF{$nc} = $avgdiff;
}

my @nc_sort_by_avgdiff = sort {$avgDIFF{$a} <=> $avgDIFF{$b}} keys %avgDIFF;

my @ref;
for (my $i=0;$i<$ctr_num;$i++){
    my $nc_file = $nc_sort_by_avgdiff[$i];
    push @ref, $nc_file;
}

# print selected refs
my $selected_ref = "$outdir/selected_ref.txt";
open O, ">$selected_ref" or die;
print O "sample\tref\tavg_diff\n";
for my $s (@ref){
    my $val = $avgDIFF{$s};
    print O "$file\t$s\t$val\n";
}
close O;

# ref matrix
my $ref_matrix_file = "$outdir/ref_matrix.xls";
open O, ">$ref_matrix_file" or die;
print O "CHR:START-END";
for (@ref){
    print O "\t$_";
}
print O "\n"; # print header: target/ref1/ref2/ref3/...

my %ref_norm_depth;
for my $s (@ref){
    open IN, "$s" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        my $target = "$arr[0]\:$arr[1]\-$arr[2]"; # chr:start-end
        push @{$ref_norm_depth{$target}}, $arr[-1];
    }
    close IN;
}

my $target_aref = &get_target($file);
for my $t (@{$target_aref}){
    my @val;
    if (exists $ref_norm_depth{$t}){
        @val = @{$ref_norm_depth{$t}};
    }else{
        die "can not find target $t in ref samples\n";
    }
    my $val = join "\t", @val;
    print O "$t\t$val\n";
}
close O;



# cal ratio/z-score
my %info;
my @ratio; # cal sampleQC

open IN, "$file" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $target = "$arr[0]\:$arr[1]\-$arr[2]"; # chr:start-end
    my $norm_depth = $arr[-1];
    my @ref_norm_depth = sort {$a <=> $b} @{$ref_norm_depth{$target}}; # sort by 1->N
    #print "ref norm depth is:\n";
    #print "@ref_norm_depth\n";

    # 去除ref首尾各5%值
    my $n_val = scalar(@ref_norm_depth);
    my $lower = floor($n_val * 0.05);
    my $upper = ceil($n_val * 0.95);

    if ($upper == $n_val){
        $upper -= 1;
    }

    my @eff_norm_depth = @ref_norm_depth[$lower..$upper];
   
    # 计算ref中位置
    my $ref_norm = &median(\@eff_norm_depth);
    my $ratio;
    if ($ref_norm != 0){
        $ratio = sprintf "%.2f", $arr[-1]/$ref_norm;
    }else{
        $ratio = "NA"; # 无效ratio
    }

    push @ratio, $ratio;

    my $zscore = &cal_zscore($arr[-1], \@eff_norm_depth);
    my $targetCV = &targetCV(\@eff_norm_depth);
    my $pval = "NA"; # 计算ref中分布的正态性
    $info{$target} = "$ratio\t$zscore\t$targetCV\t$ref_norm\t$pval"; # ratio/zscore/targetCV/ref_norm/p_val
}
close IN;


my $cnv_file = "$outdir/$name\.ratio.xls";
open O, ">$cnv_file" or die;
print O "CHR\tSTART\tEND\tGENE\tREGION_COV\tMEAN_COV\tNORMALIZED_COV\tREF_NORMALIZED_COV\tRATIO\tZSCORE\ttargetCV\tsampleCV\trefNormality\n";

# cal sampleCV
my $sampleCV = &sampleCV(\@ratio);

open IN, "$file" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $target = "$arr[0]\:$arr[1]\-$arr[2]"; # chr:start-end
    my $val = $info{$target};
    my @val = split /\t/, $val; # ratio/zscore/targetCV/ref_norm/p_val
    print O "$_\t$val[-2]\t$val[0]\t$val[1]\t$val[2]\t$sampleCV\t$val[-1]\n";
}
close IN;
close O;



####### SUB FUNC #######
sub cal_cov_diff{
    my ($file1,$file2) = @_;

    # check if file1's chr:start-end is same as file2's
    my $pos1_aref = &get_pos($file1);
    my $pos2_aref = &get_pos($file2);

    my $n1 = scalar(@{$pos1_aref});
    my $n2 = scalar(@{$pos2_aref});

    if ($n1 != $n2){
        die "the row number is not equal between $file1 and $file2\n";
    }


    my $diff_flag = 0;
    for my $i (1..$n1){
        my $idx = $i - 1;
        my $t1 = $pos1_aref->[$idx];
        my $t2 = $pos2_aref->[$idx];
        if ($t1 ne $t2){
            $diff_flag += 1;
        }
    }

    if ($diff_flag >= 1){
        # at least one target is not same
        die "at least one target(chr:start-end) is not same between $file1 and $file2\n";
    }

    my $file1_target2cov_href = &target2cov($file1);
    my $file2_target2cov_href = &target2cov($file2);

    my @abs_diff;
    for my $t (@{$pos1_aref}){
        my $val_1 = $file1_target2cov_href->{$t};
        # print "$val_1\n";
        my $val_2 = $file2_target2cov_href->{$t};
        my $abs_diff = abs($val_1 - $val_2);
        push @abs_diff, $abs_diff;
    }

    # remove 5% highest value (biggest abs diff)
    my @sort_cov_diff = sort {$b <=> $a} @abs_diff;
    my $eleNum = scalar(@abs_diff);
    my $numRemove = int(($eleNum/100)*5);
    splice(@sort_cov_diff,0,$numRemove);
    my $sum_diff = sum(@sort_cov_diff);
    my $avg_diff = sprintf "%.4f", $sum_diff/scalar(@sort_cov_diff);

    return($avg_diff);
}

sub get_target{
    my ($file) = @_;
    my @target;
    open IN, "$file" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        my $target = "$arr[0]\:$arr[1]\-$arr[2]"; # chr:start-end
        push @target, $target;
    }
    close IN;

    return(\@target);
}

sub median{
    my ($val_aref) = @_;
    my $median;
    my @val = sort {$a <=> $b} @{$val_aref};
    if (scalar(@val) % 2 == 0){
        my $v1 = $val[scalar(@val)/2-1];
        my $v2 = $val[scalar(@val)/2];
        $median = sprintf "%.2f", ($v1+$v2)/2;
    }else{
        $median = sprintf "%.2f", $val[ceil(scalar(@val)/2)];
    }

    return($median);
}

sub cal_zscore{
    my ($sampleValue, $aref) = @_;
    my ($mean,$sd) = &calMeanSD($aref);
    my $z;
    if ($sampleValue ne "NA"){
        if ($sd != 0){
            $z = sprintf "%.2f", ($sampleValue - $mean) / $sd; # 有些target在所有ref中深度均为0
        }else{
            $z = "NA";
        }
    }else{
        $z = "NA";
    }

    return($z);
}

sub targetCV{
    my ($aref) = @_;
    my ($mean,$sd) = &calMeanSD($aref);
    my $vc;
    if ($mean != 0){
        $vc = sprintf "%.2f", $sd/$mean;
    }else{
        $vc = "NA";
    }

    return($vc);

}

sub sampleCV{
    my ($aref_ratio) = @_;

    my @eff_ratio; # 去除NA
    for my $r (@{$aref_ratio}){
        next if ($r eq "NA");
        push @eff_ratio, $r;
    }

    # 去除首尾各2.5%
    my @sortedRatio = sort {$a <=> $b} @eff_ratio;
    my $numValues = scalar(@sortedRatio);
    my $low25perc = ceil((0.025*$numValues));
    my $high25perc = floor((0.975*$numValues));
    my @sliceNormVal = @sortedRatio[($low25perc) .. ($high25perc-1)];
    my ($mean,$sd) = &calMeanSD(\@sliceNormVal);
    my $sampleCV = sprintf "%.2f", $sd/$mean;

    return($sampleCV);
}

sub get_pos{ # normalized file
	my ($file) = @_;
	my @pos;
	open IN, "$file" or die;
	<IN>;
	while (<IN>){
		chomp;
		my @arr = split /\t/;
		my $region = "$arr[0]\:$arr[1]\-$arr[2]"; # chr:start-end
		push @pos, $region;
	}
	close IN;

	return(\@pos);
}
		
sub target2cov{ # normalized file
	my ($file) = @_;
	my %target2cov;
	open IN, "$file" or die;
	<IN>;
	while (<IN>){
		chomp;
		my @arr = split /\t/;
		my $target = "$arr[0]\:$arr[1]\-$arr[2]";
		$target2cov{$target} = $arr[-1];
	}
	close IN;

	return(\%target2cov);
}
	
sub calMeanSD{
	my ($values) = @_;
	my $total=0;
	my $countEle = 0;
	foreach my $ele (@$values){
		$total = $ele + $total;
		$countEle++;
	}

	my $mean = sprintf "%.4f", $total/$countEle;

	my @sqdiff;
	foreach my $ele (@$values){
		my $diff= ($ele-$mean);
		my $square= (abs($diff)*abs($diff));
		push(@sqdiff, $square);
	}

	my $totalSqDiff=0;
	foreach my $ele (@sqdiff){
		$totalSqDiff = $ele + $totalSqDiff;
	}
	my $variance = ($totalSqDiff/$countEle);
	my $sd = sqrt($variance);
	return($mean, $sd);
}







