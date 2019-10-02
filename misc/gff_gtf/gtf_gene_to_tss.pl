use strict; use warnings;

# assume one and only one gene line for each gene
# get tss from 'gene' lines in gtf, works for gencode_v29.human

while(<>){
    next if /^#/;
    my($chr, $tag, $type, $start, $end, $dot, $strand, $dot2, $anno) = split "\t", $_;
    next unless $type =~ /gene/i;

    $anno =~ /gene_id ([^;]+);/;
    my$gene_id = $1;
    $anno =~ /gene_name ([^;]+);/;
    my$gene_name = $1;

    if ($strand eq "+") {
        print "$chr $start  $start  $gene_name;$gene_id  .  $strand\n"
    }elsif ($strand eq "-") {
        print "$chr $end  $end  $gene_name;$gene_id  .  $strand\n"
    }else {warn "stand error: $_\n"}
}
