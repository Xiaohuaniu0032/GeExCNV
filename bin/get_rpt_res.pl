use strict;
use warnings;
use utf8;
use File::Basename;
use FindBin qw/$Bin/;
binmode(STDIN, ':encoding(utf8)');
binmode(STDOUT, ':encoding(utf8)');
binmode(STDERR, ':encoding(utf8)');

my ($cnv_res,$rpt_list) = @ARGV;

my $name = (split /\./, basename $cnv_res)[0];

# get exon level gene list
my $bin_dir = dirname $Bin;
my $ex_file = "$bin_dir/Ex.gene.list";

my %gene;
open IN, "$rpt_list" or die;
while (<IN>){
    chomp;
	my @arr = split /\,/, $_;
    $gene{$arr[0]} = 1;
}
close IN;

print "Gene\tTranscript\tExon\tZygosity\tClin.Signif\tSampleCV\tSampleQC_check\tCNV_quality\n";

my $qc_dir = dirname $cnv_res;
my $qc_file = "$qc_dir/$name\.ratio.xls";

die "can not find ratio file for QC" if (!-e $qc_file);

# sample qc info
my ($sampleCV,$sampleQC_check);
open IN, "$qc_file" or die;
<IN>;
my $line = <IN>;
close IN;
$sampleCV = (split /\t/, $line)[-2];

if ($sampleCV <= 0.15){
    $sampleQC_check = "PASS";
}else{
    $sampleQC_check = "Fail";
}

# cnv quality
my %exon_quality;
my $annot_file = "$qc_dir/$name\.annot.xls";
open IN, "$annot_file" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    next if (!exists $gene{$arr[4]});
    # next if ($arr[4] !~ /BRCA/);
    next if ($arr[-1] eq "NA");
    $exon_quality{$arr[4]}{$arr[-1]} = $arr[-3]; # gene=>{exon=>"cnv_quality"}
}
close IN;

my $flag = 0;
open IN, "$cnv_res" or die;
<IN>;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    if ($arr[1] =~ /BRCA/){
        $flag += 1;
        my $exon = $arr[-1];
        my $cn = $arr[-4];
        my $transcript = $arr[-2];
        my $zygosity = &zygosity($cn);
        $exon =~ s/^Exon//;
        my @exon;
        if ($exon =~ /\-/){
            # multi-exon
            my @val = split /\-/, $exon;
            for my $e ($val[0]..$val[1]){
                push @exon, $e;
            }
        }else{
            push @exon, $exon;
        }

        my @exon_quality;
        for my $e (@exon){
            my $new_ex = "exon".$e;
            my $exon_q = $exon_quality{$arr[1]}{$new_ex};
            push @exon_quality, $exon_q;
        }

        my $ex_q = join ";", @exon_quality;
        my $cnv_state;
        if ($arr[-3] =~ /del/){
            $cnv_state = "Del";
        }else{
            $cnv_state = "Amp";
        }

        my $ex_del = $arr[-1]." $cnv_state";
        print "$arr[1]\t$transcript\t$ex_del\t$zygosity\tNA\t$sampleCV\t$sampleQC_check\t$ex_q\n";
    }
}
close IN;

if ($flag == 0){
    # no brca1/2
    print "Negative\tNA\t未检测到大片段缺失/重复\tNA\tNA\t$sampleCV\t$sampleQC_check\tNA\n";
}

sub zygosity{
    my ($cn) = @_;
    my $z;
    if ($cn <= 1.3){
        $z = "Het";
    }elsif ($cn <= 0.5){
        $z = "Hom";
    }elsif ($cn >= 2.7){
        $z = "Het";
    }else{
        $z = "NA";
    }

    return($z);
}
