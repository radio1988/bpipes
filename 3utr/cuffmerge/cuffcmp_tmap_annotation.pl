#!/usr/bin/perl
#Rui LI
#Since Dr. Zhihua Jiang Lab
use strict; use warnings;

die "usage: <cuffcmp_tmap_annotation.pl> <ref_vs_query.tmap>\n" unless @ARGV == 1;
# cuffcompare -G -r query.terse.gtf -o ENS  ENS.gtf (so get all code anno for each query)
# cut -f 1,3,4,5  $PREFIX.$SAM.PAS.infos.N.rev.gtf.tmap >  $PREFIX.$SAM.PAS.infos.N.rev.gtf.tmap.txt

print "version: 2016-02-10 (use unique annotation ids)\n";


my%id2exp;#id2expression
my%id2id;
open(IN,$ARGV[0]) or die;
open(OUT,">$ARGV[0].txt") or die;

my$i = 0;
while(<IN>){
	chomp;
	$i ++;
	if($i ==1) {next}	#header
	
	my($gene_id,$loci_id,$class_code,$cuff_gene_id) = split "\t",$_;
	my$anno = $cuff_gene_id."(".$class_code.");";
	$id2exp{$gene_id} ++;
	$id2id{$gene_id} .= $anno;

}

foreach (sort keys %id2id){
	my$ids = $id2id{$_};
	my@id = split ";", $ids;
	my%uniqID;
	foreach (@id) {
		my($key,$class_code) = split (/\(/, $_);
		$class_code =~ s/\)//;
		$uniqID{$key} .= $class_code;
	}

	#get uniq id
	my$uniq_id_count; my$uniq_ids;
	foreach (keys %uniqID) {
		$uniq_id_count ++;
		$uniq_ids .= "$_ ($uniqID{$_});";
	}
		
	#output
	#print OUT "$_\t$id2exp{$_}\t$id2id{$_}\t$uniq_id_count\t$uniq_ids\n";
	print OUT "$_\t$uniq_id_count\t$uniq_ids\n";
}
