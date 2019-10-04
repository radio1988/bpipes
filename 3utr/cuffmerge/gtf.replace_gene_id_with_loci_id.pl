use strict; use warnings;

die "cmd: perl gtf.replace_gene_id_with_loci_id.pl cuffmerge.gffread_M.gtf\n" unless @ARGV==1;

while(<>){
    if (/^#/){print; next}

    # remove gene_name/gene_id, which is from fusion_loci and wrong
    $_ =~ s/gene_id[^;]+;//;
    $_ =~ s/gene_name[^;]+;//;

    # replace gene_id with loci_id
    $_ =~ s/locus "/gene_id "/;

    # remove extra spaces
    $_ =~ s/ +/ /g;

    print
}
