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
        my($gene_id) = $attribute =~ /gene_id "([^"]+)"/ || $attribute =~ /gene_id ([^;]+);/;
        my($gene_name) = $attribute =~ /gene_name "([^"]+)"/ || $attribute =~ /gene_name ([^;]+);/;
        if (!$gene_id || $gene_id eq '-'){warn "NO GENE_ID: $gene_id; $_"}
        if (!$gene_name || $gene_name eq '-'){warn "NO Gene_name: $gene_name; $_"}
        $id2name{$gene_id} = $gene_name;
    }


    if ($attribute =~ /transcript_id/){
        my($trans_id) = $attribute =~ /transcript_id "([^"]+)"/ || $attribute =~ /transcript_id ([^;]+);/;
        if (!$trans_id || $trans_id eq '-'){warn "NO trans_id: $trans_id; $_"}
    }

}

foreach my$id (sort keys %id2name){
    print "$id\t$id2name{$id}\n";
}
