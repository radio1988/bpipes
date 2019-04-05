#!/usr/bin/perl
# exact match of header in selection
# example fasta: standard, >chr1\n
# example list: >chr1\n

use strict;
use warnings;

die ('usage: <grep_fasta.pl> <xx.fa> <xxx.header.list>') unless @ARGV==2;

open(DB, $ARGV[1]) or die "err reading list $ARGV[1]\n";
open(FA, $ARGV[0]) or die "err reading fasta $ARGV[0]\n";
open(OUT, '>out.fasta') or die "err output\n";

my%accs;
while(<DB>){
    $accs{$_} = 1;
}
my$acc_size = keys %accs;
print "there are ", $acc_size, " accessions in list\n";

my$switch = 0;
while(<FA>){
    if (/^>/){
        $switch = 0;
        if ($accs{$_}){$switch = 1}
    }

    if ($switch>0){
        print OUT;
    }
}
