use strict; use warnings;

my$i=0;
while(<>){
	$_ =~ s/\r/\n/g;
	$i ++;

	if ($i == 1){
		next;
	}

	my($id, $mean, $lfc, $lfcSE, $pvalue, $fdr) = split (",", $_);
	next if $fdr eq "NA";
	next if $lfc eq "NA";

	if(($fdr <= 0.01)) {
		print "$id\n";
	}
}
