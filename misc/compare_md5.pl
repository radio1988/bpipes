#!/usr/bin/perl

use strict; use warnings;

die 'usage: perl compare_md5.pl query.md5 ref.md5' unless @ARGV == 2;

open(A, "$ARGV[0]") or die;
open(B, "$ARGV[1]") or die;

my(%query,%ref);

while(<A>){
   chmod;
   my($hash, $file) = split " ",$_;
   $query{$hash} = $file;
}

while(<B>){
    chmod;
    my($hash) = split " ", $_;
    $ref{$hash} = 1;
}

my$corrct=0;
my$wrong=0;
foreach (sort keys %query){
    if (exists $ref{$_}) {
        print join "\t", $_, $query{$_}, "validated\n";
        $corrct++;
    }
    else{
        warn join "\t", $_, $query{$_}, "not found in $ARGV[1]\n";
        $wrong++;
    }
}

print "Finished Checking, $corrct correct md5sum, $wrong wrong md5sum\n";
