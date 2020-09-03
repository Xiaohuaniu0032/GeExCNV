use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Data::Dumper;

my ($bed,$rpt_list,$geneNM,$bedtools_bin,$panel,$outdir) = @ARGV;

GetOptions(
    "bed:s" => \$bed,
    "rpt:s" => \$rpt_list,
    "generef:s" => \$geneNM,
    "bedtools:s" => \$bedtools_bin,
    "panel:s" => \$panel,
    "outdir:s" => \$outdir
    ) or die "unknown args\n";


my $cmd = "$bedtools_bin intersect -a $bed -b $geneNM -wo >$outdir/$panel\.annot.xls";
system($cmd) == 0 or die "bedtools intersect failed\n";

my %gene_NM;
# rpt_list file's format is:
# BRCA1,NM_007294
# BRCA2,NM_000059

open IN, "$rpt_list" or die;
while (<IN>){
    chomp;
    my @arr = split /\,/, $_;
    $gene_NM{$arr[0]} = $arr[1];
}
close IN;

open O, ">$outdir/$panel\.annot.rpt.xls";
my $annot = "$outdir/$panel\.annot.xls";
open IN, "$annot" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    if (exists $gene_NM{$arr[3]}){
        if ($arr[-3] eq $gene_NM{$arr[3]}){
            print O "$_\n";
        }
    }
}
close IN;
close O;



