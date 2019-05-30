# for Dr. Nathan Lawson
# Converting GTF to UCSC format
# change from 1 to chr
# KZ116046 to chrUn_KZ116046v1
# etc.
print "usage: perl gtf_chrName_change.pl GCF_000002035.6_GRCz11_genomic.gff ucscToRefSeq.txt\n";
use strict; 
use warnings;

die "Usage: perl xxx.pl xxx.gtf ucscToRefSeq.txt\n" unless @ARGV==2;

open(GTF, $ARGV[0]) or die "err reading GTF\n";
open(DB, $ARGV[1]) or die "err reading chromInfo.txt\n";

my%hash;
while(<DB>){
	$_ =~ s/[\r\n]//g;
    my @ele = split('\t', $_);
    $hash{$ele[3]} = $ele[0]
}
close DB;

open (TAB,">hash.txt") or die "err outputing hash\n";
foreach my$key (sort keys %hash){
	print TAB "$key\t$hash{$key}\n";
}
close TAB;

open(OUT, ">output.gtf") or die "err output\n";
open(ERR, ">err.gtf") or die "err stderr\n";

while(<GTF>){
	if ($_ =~ m/^#/){print OUT; next;}  # header
	$_ =~ s/[\r\n]//g;

        my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $_, 9);

	if($hash{$seqname}){
		$seqname = $hash{$seqname};
        	my$out_gtf = join("\t", $seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute);
		print OUT "$out_gtf\n";
	}else{
        	my$out_gtf = join("\t", $seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute);
		print ERR "$out_gtf\n";
	}
}
close GTF;
close OUT;
close ERR;
