#!/usr/bin/perl
use strict; use warnings;

die ("usage: <perl> <vcf.get_mut_vs_wt.pl> <xx.snpEff.cancer.vcf>\n") unless @ARGV==1;
print "get mutations that are in the mut, but not in wt from the VCF
the input VCF should be the output of snpEff cancer pipeline\n";

my$in = $ARGV[0];
my$out = $in;
$out =~ s/vcf$/mut_vs_wt\.vcf/;

open(IN, $ARGV[0]) or die;
open(OUT, ">$out") or die;

my$count=0;
my$lines=0;
while (<IN>){
    my$line = $_;

    if ($line =~ m/^#/){print OUT}

    $lines++;

    if ($line =~ m/ANN=[^\|]*-/){
        print OUT $line;
        $count++;
    }
}

print "there are $count mut_vs_wt among $lines of input $in\n";
