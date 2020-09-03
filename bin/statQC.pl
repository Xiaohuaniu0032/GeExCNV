use strict;
use warnings;
use POSIX qw(ceil);
use POSIX qw(floor);
use List::Util qw(sum);
# 统计样本QC

my ($bam,$name,$norm_depth,$samtools_bin,$outdir) = @ARGV;

my $outfile = "$outdir/$name\.QC.xls";
open O, ">$outfile" or die;
print O "name\treads_num\tproper_aln_pct\tmean_cov\t50X\t100X\t200X\t300X\t400X\t500X\t0.2X\t0.5X\tsampleQC\n";


# 统计flagstat信息
my $flagstat_file = "$outdir/$name\.flagstat.txt";
print "$flagstat_file\n";

my $cmd = "$samtools_bin flagstat $bam >$flagstat_file";

system($cmd) == 0 or die "samtools flagstat failed\n";

my ($reads_n,$aln_pct);
open IN, "$flagstat_file" or die "can not find $flagstat_file\n";
while (<IN>){
    chomp;
    if (/paired in sequencing/){
        $reads_n = (split /\s/, $_)[0];
    }
    if (/properly paired/){
        my $aln_num = (split /\s/, $_)[0];
        $aln_pct = sprintf "%.2f", $aln_num/$reads_n;
    }
}
close IN;

my ($mean_cov,$x02,$x05,$x50,$x100,$x200,$x300,$x400,$x500,$sample_qc);
my @norm_depth;
my @depth;
open IN, "$norm_depth" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    $mean_cov = $arr[-2];
    push @norm_depth, $arr[-1];
    push @depth, $arr[-3];
}
close IN;

$x02 = &covX(\@depth,int(0.2*$mean_cov));
$x05 = &covX(\@depth,int(0.5*$mean_cov));
$x50 = &covX(\@depth,50);
$x100 = &covX(\@depth,100);
$x200 = &covX(\@depth,200);
$x300 = &covX(\@depth,300);
$x400 = &covX(\@depth,400);
$x500 = &covX(\@depth,500);

my @raitoF = glob "$outdir/*.ratio.xls";
my $raitoF = $raitoF[0];
open IN, "$raitoF" or die;
<IN>;
my $line = <IN>;
close IN;
my $sampleqc = (split /\t/, $line)[-2];

print O "$name\t$reads_n\t$aln_pct\t$mean_cov\t$x50\t$x100\t$x200\t$x300\t$x400\t$x500\t$x02\t$x05\t$sampleqc\n";
close O;



sub covX{
    my ($cov_aref,$cutoff) = @_;
    my $num = 0;
    my $n;
    for my $d (@{$cov_aref}){
        $n += 1;
        if ($d >= $cutoff){
            $num += 1;
        }
    }
    my $pct = sprintf "%.2f", $num/$n;
    return($pct);
}





