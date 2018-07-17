use strict; use warnings;

my%id2name;
while(<>){
    if (/^#/){next}

    $_ =~ s/[\r\n]//g;

    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $_, 9);


    if ($attribute =~ /gene_name/ && $attribute =~ /gene_id/){
        my($gene_id) = $attribute =~ /gene_id "([^"]+)"/;
        if (!$gene_id){warn $_}
        my($gene_name) = $attribute =~ /gene_name "([^"]+)"/;
        $id2name{$gene_id} = $gene_name;
    }

}

foreach my$id (sort keys %id2name){
    print "$id\t$id2name{$id}\n";
}
