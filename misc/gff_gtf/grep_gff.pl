#!/usr/bin/perl
# exact match of header in selection
# example gff: stanard, chr1, chr10
# example list: 
#   >chr1
#   >chr10
#   >chr11

use strict;
use warnings;

die ('usage: <grep_fasta.pl> <xx.gff/gtf> <xxx.fasta.header>') unless @ARGV==2;

open(DB, $ARGV[1]) or die "err reading list $ARGV[1]\n";
open(GFF, $ARGV[0]) or die "err reading gff/gtf $ARGV[0]\n";
open(OUT, '>out.gff') or die "err output\n";

my%accs;
while(<DB>){
    chomp;
    $_ =~ s/^>\s*//;
    $accs{$_} = 1;
}
my$acc_size = keys %accs;
print "there are ", $acc_size, " accessions in list\n";

my$switch = 0;
while(<GFF>){
    if (/^#/){print OUT; next}

    my($chr, @others) = split "\t", $_;

    if ($accs{$chr}){print OUT}
}
