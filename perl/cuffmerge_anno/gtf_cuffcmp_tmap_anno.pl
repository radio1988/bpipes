use strict; 
use warnings;

die "Usage: <gtf_cuffcmp_tmap_anno.pl> <input.gtf> <vs_ref.tmap> <ref.gtf>" unless @ARGV==3;

open(GTF, $ARGV[0]) or die 'err reading GTF';
open(TMAP, $ARGV[1]) or die "err reading tmap\n";
open(GTFREF, $ARGV[2]) or die "err reading GTFREF\n";

## Build Lookup Tables
my %tr2genename;
my %tr2transid; 
my %tr2classcode;
while(<TMAP>){
	chomp;
    if (/^ref_gene_id/){next}  #header
	my($ref_gene_name, $ref_tr_id, $class_code, $in_gene_id, $in_tr_id) = split "\t", $_;
	$tr2genename{$in_tr_id} = $ref_gene_name;
	$tr2transid{$in_tr_id} = $ref_tr_id; 
	$tr2classcode{$in_tr_id} = $class_code;
}
close TMAP;

## Tr_ID to Gene_ID
my%reftr2regene;
while(<GTFREF>){
    if (/^#/){next}

    $_ =~ s/[\r\n]//g;
   
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $_, 9);

    if ($attribute !~ /transcript_id/){
        # some lines don't have trans_id
        next;
    }
    
    my($tr_id) = $attribute =~ /transcript_id "([^"]+)"/;
    my($gene_id) = $attribute =~ /gene_id "([^"]+)"/;
    $reftr2regene{$tr_id} = $gene_id;
}
close GTFREF;

## Output Anno Info
open(TAB, ">$ARGV[1].txt") or die "err output lookup table\n";
print TAB "in_trans_id ref_trans_id ref_gene_id ref_gene_name    tmap_class_code\n";
foreach my$in_tr(sort keys %tr2genename){
    my$ref_trid = $tr2transid{$in_tr};
    my$ref_geneid = $reftr2regene{$ref_trid} || '-';
    my$ref_gene_name = $tr2genename{$in_tr};
    my$class_code = $tr2classcode{$in_tr};
	print TAB ("$in_tr  $ref_trid   $ref_geneid $ref_gene_name  $class_code\n");
}
close TAB;

## Annotate and Transform GTF
open(OUT, ">$ARGV[0].$ARGV[1].annotated") or die "err output anno GTF\n";
open(ERR, ">$ARGV[0].$ARGV[1].nonanno") or die "err output to nonanno GTF\n";
while(<GTF>){
    if ($_ =~ m/^#/){print OUT; next;}  # header

    $_ =~ s/[\r\n]//g;

    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $_, 9);
    my($tr_id) = $attribute =~ /transcript_id "([^"]+)"/;	
    my($exon) = $attribute =~ /exon_number "([^"]+)"/;	

    my$ref_trans_id = $tr2transid{$tr_id};
    my$gene_id = $reftr2regene{$ref_trans_id};

    my$genename = $tr2genename{$tr_id};
    my$classcode = $tr2classcode{$tr_id};
    if ($classcode =~ /[=cjeo]/){
        my$attribute_new = "gene_id \"$gene_id\"; transcript_id \"$tr_id\"; exon_number \"$exon\"; gene_name \"$genename\"; class_code \"$classcode\";";
        my$out = join("\t", $seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute_new);	
        print OUT $out, "\n";
    }else{
        my$attribute_new = $attribute;
        my$out = join("\t", $seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute_new);
        print ERR $out, "\n";
    }
}
close GTF;
close OUT;
close ERR;
