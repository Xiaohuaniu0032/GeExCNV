use strict;
use warnings;

# http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1220235277_SJWvNhqI4wCU4FbY6LcxOLoaaC5l&hgta_doSchemaDb=macFas5&hgta_doSchemaTable=refFlat
# Note: all start coordinates in our database are 0-based, not 1-based. See explanation here.

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

	# add cds region
	my $cds_start = $arr[6] + 1;
	my $cds_end   = $arr[7];

	if ($strand eq "+"){
		for (my $i=1;$i<=$exon_num;$i++){
			my $start = $start_exon[$i-1] + 1; # 0-based
			my $end = $end_exon[$i-1];
			$exon_idx += 1;
			my $exon = "exon".$exon_idx;
			print O "$arr[2]\t$start\t$end\t$arr[0]\t$strand\t$arr[1]\t$exon\tCDS:$cds_start\:$cds_end\n";
		}
	}else{
		for (my $i=1;$i<=$exon_num;$i++){
			my $start = $start_exon[$i-1] + 1; # last exon -> first exon
			my $end = $end_exon[$i-1];
			$exon_idx = $exon_num - $i + 1;
			my $exon = "exon".$exon_idx;
			print O "$arr[2]\t$start\t$end\t$arr[0]\t$strand\t$arr[1]\t$exon\tCDS\:$cds_start\:$cds_end\n";
		}
	}
}
close IN;
close O;
	

	
