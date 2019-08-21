use strict; use warnings;

print "gtf_2_anno_tab.gene_level.pl <xxx.gtf>\n";
print "Out: anno.txt\n";
print "good for ENSEMBL MOUSE v94 ANNOTATION\n";

die unless @ARGV == 1;

open(GTF, shift) or die;

my%id2name;
while(<GTF>){
    chomp;
    next if /^#/;
    next if $_ !~ m/gene_name/;
    next if $_ !~ m/gene_id/;

    my($name) = $_ =~ /gene_name([^;]+);/;
    my($id) = $_ =~ /gene_id ([^;]+);/;
    my($type) = $_ =~ /gene_biotype ([^;]+);/;
    if (length($type) == 0){$type = "-"}

    $id2name{$id} = "$name\t$type";
}
close GTF;

open(OUT, ">anno.txt");
print OUT "Gene\tName\tType\n";
for (sort keys %id2name){
    print OUT "$_\t$id2name{$_}\n"
}
close OUT;
