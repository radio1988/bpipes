use strict; use warnings;

my$inname = $ARGV[0];
my$inprefix = $inname;
$inprefix =~ s/txt//;

print("For $inname\n");

open(IN, $inname) or die;

my$nblocks=0;
while (<IN>){
	if (/ for My Projects->Amelia->/){
		my$nblocks ++;
		if ($nblocks > 1) {close OUT}

		my($a, $b) = /(^.*) for My Projects->Amelia->(.*)$/;
		$a =~ s/ /_/g;
		
		my$outname = $inprefix.$a.".s.txt";
		open(OUT, ">$outname") or die;
		print OUT "# $_";
	} else {
		print OUT;
	}
}

print("\n\n");
