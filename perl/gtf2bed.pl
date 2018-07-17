#!/usr/bin/perl
# from julie.v14 to bed for lift-over  
use strict;
use warnings;

while(<>){
	$_ =~ s/[\r\n]//g;
	my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $_, 9);
	my ($trans_id) = $attribute =~ m/(transcript_id [^;]+;)/;
	$attribute =~ s/\s/%/g;  # keep attri
	## output gtf
	my $bed = join("\t", $seqname, $start-1, $end-1, ($source.';'.$feature.';'.$attribute), $score, $strand, );
	print "$bed\n";

}
