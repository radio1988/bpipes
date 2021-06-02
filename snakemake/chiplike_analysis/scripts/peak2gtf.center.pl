#!/bin/perl
# input narrowPeak/bed/gappedPeak
# output gtf for featureCount

while(<>){
    next if /^#/;
    chomp;
    my@F=split "\t", $_;
    my$chr=$F[0];
    my$start=$F[6];
    my$end=$F[7];
    my$peak_id=$F[3];
    my$score = ".";
    my$strand = ".";
    print join "\t", $chr, "PeakCaller", "peak", $start+1, $end+1, ".", ".", ".", "gene_id \"$peak_id\";\n"

}
