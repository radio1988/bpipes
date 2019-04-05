#!/usr/bin/perl
use strict;
use warnings;
# from bed to gtf after lift-over, format paired with gtf2bed(special attri)
while(<>){
	$_ =~ s/[\r\n]//g;
	my ($chr, $start, $end, $name, $score, $strand, $extras) = split(/\t/, $_, 7);
	## replace $chr
	if ($chr =~ m/^chr\d+$/){
		# chr4 -> 4
		$chr =~ s/^chr//g;
	}
	elsif ($chr =~ m/^chrUn_.+v1$/){
		# chrUn_KN150269v1 -> KN150269
		$chr =~ s/^chrUn_//g;
		$chr =~ s/v1$//g;
	}

	## parse name (nathen v14 zebrafish specific)
	my ($source, $feature, $attribute) = split(';', $name, 3);
	$attribute =~ s/%/ /g;

	## output gtf
	my $gtf = join("\t", $chr, $source, $feature, $start+1, $end+1, $score, $strand, '.', $attribute );
	print "$gtf\n";

}
