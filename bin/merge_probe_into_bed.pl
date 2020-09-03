use strict;
use warnings;

my ($in) = @ARGV;
my %info;
open IN, "$in" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $val = "$arr[1]\t$arr[2]";
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
        push @{$info{$chr_int}}, $val;
    }else{
        push @{$info{$arr[0]}}, $val;
    }
}
close IN;

my %bed;
foreach my $c (sort {$a <=> $b} keys %info){
    my $chr;
    if ($c == 23){
        $chr = "chrX";
    }elsif ($c == 24){
        $chr = "chrY";
    }else{
        $chr = "chr".$c;
    }

    my @val = @{$info{$c}};
    my @seg;
    my $first = shift @val;
    push @seg, $first;
    for my $val (@val){
        my $before_val = $seg[-1];
        my @before_val = split /\t/, $before_val;

        my @this_val = split /\t/, $val;
        if ($this_val[0] <= $before_val[-1]){
            my $new_seg = "$before_val[0]\t$this_val[-1]";
            pop @seg;
            push @seg, $new_seg;
        }else{
            push @seg, $val;
        }
    }

    for my $s (@seg){
        my @val = split /\t/, $s;
        print "$chr\t$val[0]\t$val[1]\n";
    }
}



