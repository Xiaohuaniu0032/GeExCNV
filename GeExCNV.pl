use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw/$Bin/;
use List::Util qw(sum);
use POSIX qw(ceil);
use POSIX qw(floor);

# 2019-2020

my ($bam,$cfg,$rpt_gene,$outdir);
GetOptions(
    "bam:s" => \$bam,          # NEED
    "cfg:s" => \$cfg,          # NEED
    "rpt:s" => \$rpt_gene,     # NEED
    "outdir:s" => \$outdir,    # NEED
    ) or die "unknown args\n";


if (not defined $bam || not defined $cfg || not defined $outdir){
    die "please check your args\n";
}

# check files
die "can not find $bam file" if (!-e $bam);
die "can not find $cfg file" if (!-e $cfg);

if (!-d $outdir){
    `mkdir -p $outdir`;
}

# parse cfg.txt and print cfg info into a cfg.log file
my %cfg;

my $cfg_log = "$outdir/cfg.log";
open CFG, ">$cfg_log" or die;
print CFG "####### config info #######\n";

open IN, "$cfg" or die;
while (<IN>){
    chomp;
    next if /^$/;
    next if (/^\#/);
    my @val = split /\:/;
    print CFG "$val[0] ===> $val[1]\n";
    $cfg{$val[0]} = $val[1];
}
close IN;
close CFG;


# auto detect sample name
my $name = (split /\./, basename $bam)[0];
print "auto detected sample name is: $name\n";

my $runsh = "$outdir/cnv\_$name\.sh";
open O, ">$runsh" or die;
print O "cd $outdir\n";

# s1. cal depth and mean norm
print O "\# step1:cal depth\n";
my $bed = $cfg{"bed"};
my $samtools = $cfg{"samtools"};
print O "perl $Bin/bin/cal_depth.pl -bam $bam -bed $bed -samtools $samtools -outdir $outdir\n\n";

# s2. cal copy number ratio
print O "\# step2:cal copy number ratio\n";
my $norm_depth = "$outdir/$name\.normalized_depth.txt";
my $ctr_num = $cfg{"ctr_num"};
my $ctr_dir = $cfg{"ctr_dir"};
print O "perl $Bin/bin/cal_cnr.pl -infile $norm_depth -ctr_num $ctr_num -ctr_dir $ctr_dir -outdir $outdir\n\n";

# s3. annot NM/exon and cnv state(Del/Amp/Norm)
print O "\# step3:annot exon\n";
my $ratioF = "$outdir/$name\.ratio.xls";
my $annotF = $cfg{"annot_file"};
my $annot_of = "$outdir/$name\.annot.xls";
print O "perl $Bin/bin/annot_exon_CNV_state.pl -in $ratioF -o $annot_of -annot $annotF\n\n";

# s4. get qc info
print O "\# step4:cal qc\n";
print O "perl $Bin/bin/statQC.pl $bam $name $norm_depth $samtools $outdir\n\n";

# s5. get exon-level and gene-level CNV result
my $strandF = $cfg{"strand_file"};
print O "perl $Bin/bin/mergeCNV.pl -n $name -in $annot_of -bed $bed -annot $annotF -strand $strandF -rpt $rpt_gene\n\n";

# s6. get rpt result
my $exon_cn = "$outdir/$name\.exon.CNV.xls";
my $gene_cn = "$outdir/$name\.gene.CNV.xls";
my $result_file = "$outdir/$name\.rpt.cnv.xls";
print O "perl $Bin/bin/get_rpt_res.pl $exon_cn $rpt_gene >$result_file\n\n";

# s7. plot fig
my $fig_cn_pdf = "$outdir/$name\.CNV.pdf";
my $rscript=$cfg{"Rscript"};
print O "$rscript $Bin/bin/plot_copy_num.R $ratioF $fig_cn_pdf\n\n";
print O "/usr/bin/convert -quality 100 -density 360 $fig_cn_pdf $outdir/$name\.CNV.png\n\n";

close O;



