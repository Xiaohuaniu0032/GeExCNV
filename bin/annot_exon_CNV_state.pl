use strict;
use warnings;
use Getopt::Long;
use FindBin qw/$Bin/;

my ($infile,$del_cutoff_low,$del_cutoff_high,$amp_cutoff_low,$amp_cutoff_high,$targetQC,$sampleQC,$zscore_cutoff_low,$zscore_cutoff_high,$outfile);
my ($annot_file);

GetOptions(
	"in:s" => \$infile,                       # need
	"del_low:f" => \$del_cutoff_low,          # default: 0.6
	"del_high:f" => \$del_cutoff_high,        # default: 0.65
	"amp_low:f" => \$amp_cutoff_low,          # default: 1.35
	"amp_high:f" => \$amp_cutoff_high,        # default: 1.45
	"targetQC:f" => \$targetQC,               # default: 0.2
	"sampleQC:f" => \$sampleQC,               # default: 0.3
	"z_low:f" => \$zscore_cutoff_low,         # default: 3
	"z_high:f" => \$zscore_cutoff_high,       # default: 4
	"o:s" => \$outfile,                       # need
	"annot:s" => \$annot_file,                # need
	) or die "unkonw args\n";


# 默认值
if (not defined $del_cutoff_low){
	$del_cutoff_low = 0.6;
}
if (not defined $del_cutoff_high){
	$del_cutoff_high = 0.65;
}

if (not defined $amp_cutoff_low){
	$amp_cutoff_low = 1.35;
}
if (not defined $amp_cutoff_high){
	$amp_cutoff_high = 1.45;
}

if (not defined $zscore_cutoff_low){
	$zscore_cutoff_low = 3;
}
if (not defined $zscore_cutoff_high){
	$zscore_cutoff_high = 4;
}

if (not defined $sampleQC){
	$sampleQC = 0.15;
}

if (not defined $targetQC){
	$targetQC = 0.15;
}

if (not defined $infile or not defined $outfile or not defined $annot_file){
	die "some must given not find\n";
}


# 需要注释基因的转录本/exon号
my %target2exon;
open IN, "$annot_file" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $target = "$arr[0]\t$arr[1]\t$arr[2]"; # chr/start/end
	$target2exon{$target} = "$arr[-2]\t$arr[-3]"; # exon/transcript
}
close IN;


open IN, "$infile" or die;
my $h = <IN>;
close IN;
chomp $h;
open O, ">$outfile" or die;
print O "sample\t$h\tsampleQC_check\ttargetQC_check\tCNV_state\tCNV_Quality\tTranscript\tExon\n";

my ($sampleQC_info,$targetQC_info);
my ($cnv_state,$cnv_quality);

open IN, "$infile" or die;
<IN>;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $target = "$arr[0]\t$arr[1]\t$arr[2]"; # chr/start/end
	my @ar = split(/\//,"$infile");
	my $smp_file = $ar[-1];
	@ar = split(/\./,"$smp_file");
	my $smp = $ar[0];

	
	# 外显子号/转录本号
	my ($exon,$NM);
	if (exists $target2exon{$target}){
		($exon,$NM) = split /\t/, $target2exon{$target};
	}else{
		($exon,$NM) = ("NA","NA");
	}

	# sampleQC
	if ($arr[-2] <= $sampleQC){
		$sampleQC_info = "PASS";
	}else{
		$sampleQC_info = "Fail";
	}

	# ratio/zscore/targetCV有可能为NA
	if ($arr[-5] eq "NA" || $arr[-4] eq "NA" || $arr[-3] eq "NA"){
		$targetQC_info = "NA";
		$cnv_state = "NA";
		$cnv_quality = "NA";
	}else{
		# targetQC
		if ($arr[-3] <= $targetQC){
			$targetQC_info = "PASS";
		}else{
			$targetQC_info = "Fail";
		}

		# 先使用宽松的条件判断exon状态
		if ($arr[8] <= $del_cutoff_high and abs($arr[9]) >= $zscore_cutoff_low){
			# ratio <= 0.65 && zscore <= -3
			$cnv_state = "Del";
		}elsif($arr[8] >= $amp_cutoff_low and $arr[9] >= $zscore_cutoff_low){
			# ratio >= 1.35 && zscore >= 3
			$cnv_state = "Amp";
		}else{
			$cnv_state = "Normal";
		}

		# 判断可信度
		if ($cnv_state eq "Del"){
			# for Del
			if ($arr[8] <= $del_cutoff_low and abs($arr[9]) >= $zscore_cutoff_high){
				# ratio和zscore都高可信
				$cnv_quality = "High";
			}elsif ($arr[8] <= $del_cutoff_low and (abs($arr[9]) >= $zscore_cutoff_low and abs($arr[9]) <= $zscore_cutoff_high)){
				# ratio高可信
				$cnv_quality = "Middle";
			}elsif (($arr[8] >= $del_cutoff_low and $arr[8] <= $del_cutoff_high) and abs($arr[9]) >= $zscore_cutoff_high){
				# zscore高可信
				$cnv_quality = "Middle";
			}elsif (($arr[8] >= $del_cutoff_low and $arr[8] <= $del_cutoff_high) and (abs($arr[9]) >= $zscore_cutoff_low and abs($arr[9]) <= $zscore_cutoff_high)){
				$cnv_quality = "Low";
			}
		}else{
			# for Amp
			$cnv_quality = "NA";
		}
	}

	print O "$smp\t$_\t$sampleQC_info\t$targetQC_info\t$cnv_state\t$cnv_quality\t$NM\t$exon\n";
}

close IN;
close O;
