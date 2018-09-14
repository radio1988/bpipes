#!/usr/bin/perl
use strict; use warnings;

sub is_ref{
    if ($_[0] eq '0/0'){
        return 1
    } else {
        return 0;
    }
}

sub is_mut{
    if ($_[0] ne '0/0' && $_[0] ne './.'){
        return 1
    } else {
        return 0;
    }
}

sub contain_mut{
    foreach (@_) {
        if ($_ ne '0/0' && $_ ne './.'){
            return 1;
        }
    }
    return 0;
}

sub contain_nwt{
	my$WT= splice @_, 0, 1;
	my@MUTS = @_;
	foreach (@MUTS){
		if (($_ ne './.') && ($_ ne $WT)){
			return 1;
		}
	}
	return 0;
}


die ("usage: <perl> <vcf.get_mut.pl> <xx.samples_merged.vcf> <WT_ID>\n") unless @ARGV==2;

print "get mutations that are 0/0 in WT, but not all 0/0 in MUTs\n";
print "OR: get mutations that is_mut in WT, but not all WT in MUTs\n";
print "./. not counted as MUT\n";

print "WT_ID should be 1-#sample, the index of WT sample\n";

my$in = $ARGV[0];
open(IN, $ARGV[0]) or die;

my$out0 = $in;
$out0 =~ s/vcf$/mut0\.vcf/;
print "out0: $out0 contains variants with 0/0 in wt\n";
open(OUT0, ">$out0") or die;


my$out1 = $in;
$out1 =~ s/vcf$/mut1\.vcf/;
print "out1: $out1 contains variants with non-ref in wt\n";
open(OUT1, ">$out1") or die;


my$WT_ID=$ARGV[1]-1;


my$count=0;
my$lines=0;
while (<IN>){
    chomp;
    my$line = $_;

    if ($line =~ m/^#/){
        print OUT0 "$line\n";
        print OUT1 "$line\n";
        next
        }

    $lines++;

    # filter
    my@eles = (split "\t", $_);
    my@standard = splice @eles, 0, 9;
    my@samples = @eles;

    my$WT = splice @samples, $WT_ID, 1;
    ($WT) = $WT =~ m/(^[^:]+):/;

    my@MUTS = @samples;
    foreach my$MUT (@MUTS) {
        ($MUT) = $MUT =~ m/(^[^:]+):/;
    }

    if ( is_ref($WT) && contain_mut(@MUTS))
    {
        print OUT0 "$line\n";
        $count++;
    } elsif ( is_mut($WT) && contain_nwt( ($WT, @MUTS) ) )
    {
        print OUT1 "$line\n";
        $count++;
    }
}

print "there are $count muts among $lines of input $in\n";
