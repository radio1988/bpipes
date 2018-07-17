use strict; use warnings;

warn "extract gene_id to gene_name table for DESeq annotation";
warn "also does GTF integrety check: valid gene_id, gene_name, trans_id\n";
warn "need quotes in ID/Names\n";

my%id2name;
while(<>){
    if (/^#/){next}

    $_ =~ s/[\r\n]//g;

    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $_, 9);


    if ($attribute =~ /gene_name/ && $attribute =~ /gene_id/){
        my($gene_id) = $attribute =~ /gene_id "([^"]+)"/;
        my($gene_name) = $attribute =~ /gene_name "([^"]+)"/;
        my($trans_id) = $attribute =~ /transcript_id "([^"]+)"/;
        if (!$gene_id || $gene_id eq '-'){warn "NO GENE_ID: $_"}
        if (!$gene_name || $gene_name eq '-'){warn "NO GENE_ID: $_"}
        if (!$trans_id || $trans_id eq '-'){warn "NO GENE_ID: $_"}
        $id2name{$gene_id} = $gene_name;
    }

}

foreach my$id (sort keys %id2name){
    print "$id\t$id2name{$id}\n";
}
