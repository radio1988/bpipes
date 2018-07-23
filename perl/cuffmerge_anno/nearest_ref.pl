use strict;
use warnings;

while(<>){
	$_ =~ m/nearest_ref \"([^\"]+)\"/;
	my $nearest_ref = $1;
	print "$nearest_ref\n";
}
