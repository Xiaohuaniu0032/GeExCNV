use strict;
use warnings;
use File::Basename;

my ($dir) = @ARGV;

my @file = glob "$dir/*/*.rpt.cnv.xls";
my $f1 = $file[0];
open IN, "$f1" or die;
my $h = <IN>;
chomp $h;
close IN;
print "name\t$h\n";

for my $f (@file){
	my $name = (split /\./, basename $f)[0];
	open IN, "$f" or die;
	<IN>;
	while (<IN>){
		chomp;
		print "$name\t$_\n";
	}
	close IN;
}

