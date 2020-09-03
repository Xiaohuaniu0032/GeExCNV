 use strict;
 use warnings;

 my ($in) = @ARGV;
 my %info;
 open IN, "$in" or die;
 while (<IN>){
    chomp;
    next if /^$/;
    my @arr = split /\t/;
    if ($arr[0] =~ /chr/){
        my $chr = $arr[0];
        $chr =~ s/^chr//;
		my $chr_int;
		if ($chr eq "X"){
			$chr_int = 23;
		}elsif ($chr eq "Y"){
			$chr_int = 24;
		}else{
			$chr_int = $chr;
		}
        $info{$chr_int}{$arr[1]} = $arr[2]; # chr/start/end
    }else{
        $info{$arr[0]}{$arr[1]} = $arr[2];
    }
 }
 close IN;

 foreach my $c (sort {$a <=> $b} keys %info){
    my $chr;
	if ($c == 23){
		$chr = "chrX";
	}elsif ($c == 24){
		$chr = "chrY";
	}else{
		$chr = "chr".$c;
	}
    my @start = sort {$a <=> $b} keys %{$info{$c}};
    for my $s (@start){
        my $end_pos = $info{$c}{$s};
        print "$chr\t$s\t$end_pos\n";
    }
 }
