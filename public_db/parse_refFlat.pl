use strict;
use warnings;

# 解析refFlat文件
my ($refFlat, $outfile) = @ARGV;

open O, ">$outfile" or die;
#print O "chr\tstart\tend\tgene\tstrand\ttranscript\texon\n";

open IN, "$refFlat" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $strand = $arr[3];
	my $start_exon = $arr[-2];
	$start_exon =~ s/\,$//;
	my $end_exon = $arr[-1];
	$end_exon =~ s/\,$//;
	
	my @start_exon = split /\,/, $start_exon;
	my @end_exon = split /\,/, $end_exon;

	my $exon_num = $arr[-3];
	my $exon_idx = 0;
	if ($strand eq "+"){
		for (my $i=1;$i<=$exon_num;$i++){
			my $start = $start_exon[$i-1];
			my $end = $end_exon[$i-1];
			$exon_idx += 1;
			my $exon = "exon".$exon_idx;
			print O "$arr[2]\t$start\t$end\t$arr[0]\t$strand\t$arr[1]\t$exon\n";
		}
	}else{
		for (my $i=1;$i<=$exon_num;$i++){
			my $start = $start_exon[$i-1];
			my $end = $end_exon[$i-1];
			$exon_idx = $exon_num - $i + 1;
			my $exon = "exon".$exon_idx;
			print O "$arr[2]\t$start\t$end\t$arr[0]\t$strand\t$arr[1]\t$exon\n";
		}
	}
}
close IN;
close O;
	

	
