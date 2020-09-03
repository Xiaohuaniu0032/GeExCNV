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

my ($bam,$bed,$samtools_bin,$outdir);

GetOptions(
    "bam:s" => \$bam,                     # NEED
    "bed:s" => \$bed,                     # NEED
    "samtools:s" => \$samtools_bin,       # NEED
    "outdir:s" => \$outdir,               # NEED
    ) or die "unknown args\n";

# check args
die "can not find bam or bed file(s), please check args" if (!-e $bam || !-e $bed);

my $name = (split /\./, basename $bam)[0];
print "auto detected sample name (from bam) is: $name\n";

print "cal cov depth for $name\n";
# cal depth
&countDepth($bam,$bed,$outdir);

# mean depth normalization
my $cov_file = "$outdir/depth.txt";
my $norm_file = "$outdir/$name\.normalized_depth.txt";
&normalize_depth($cov_file,$norm_file);


####### SUB FUNC #######
sub countDepth{
    my ($bam, $bed, $outdir) = @_;
    my $depth_file = "$outdir/bedcov.txt";
    my $cmd = "$samtools_bin bedcov $bed $bam >$depth_file"; # bedcov will skip markdup reads automatically
    system($cmd) == 0 or die "samtools bedcov failed\n";

    # cal mean cov
    my $cov = "$outdir/depth.txt";
    open O, ">$cov" or die;
    print O "chr\tstart\tend\tgene\tdepth\n";

    open IN, "$depth_file" or die;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        next if ($arr[0] =~ /[xXyY]/); # skip X/Y chr
        my $len = $arr[2] - $arr[1];
        my $cov = int($arr[-1]/$len);
        print O "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$cov\n";
    }
    close IN;

    `rm $depth_file`;
}


sub normalize_depth{
    my ($infile,$outfile) = @_;

    my @cov;
    open IN, "$infile" or die;
    my $h = <IN>;
    chomp $h;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        push @cov, $arr[-1];
    }
    close IN;

    my $sum_cov;
    for (@cov){
        $sum_cov += $_;
    }

    my $mean_cov = int($sum_cov/scalar(@cov));
    
    open O, ">$outfile" or die;
    print O "$h\tmean_depth_autochr\tnormalized_depth\n";

    open IN, "$infile" or die;
    <IN>;
    while (<IN>){
        chomp;
        my @arr = split /\t/;
        my $norm_cov = sprintf "%.2f", $arr[-1]/$mean_cov;
        print O "$_\t$mean_cov\t$norm_cov\n";
    }
    close IN;
    close O;
}


