# find table from gene_id to gene_name
# find cases multiple gene_id has same gene_name
# output ID to GeneID, ID to GeneName
# perl gff_GeneID2Name.pl GCF_000002035.6_GRCz11_genomic.ucsc.primary.gff > ID2GeneID.txt
use strict; use warnings;
die 'usage: <gtf_GeneID2Name.pl> <ncbi.refseq.gff3> ' unless @ARGV == 1;

my %GeneID2ID;
my %ID2GeneID;


my$i=0;
while(<>) {
	if (/^#/){
		unless(/^##/){print}  # redundant things
		next;
	}
	if (length($_) < 1) {print; next}
	$i++;
	# die if $i > 10;

	my(@F) = split "\t", $_;
	# skipping list (without GeneID)
	if ($F[2] eq 'region'){next}
	if ($F[2] eq 'cDNA_match'){next}
	if ($F[2] eq 'match'){next}
	if ($F[2] eq 'D_loop'){next}

	# Get GeneID, no lines without GeneID
	my$GeneID = "";
	if (/GeneID:/){
		if (m/GeneID:([^,]+),/) {$GeneID = $1} 
		elsif (m/GeneID:([^"]+)"/) {$GeneID = $1}
		elsif (m/GeneID:([^;]+);/) {$GeneID = $1}
		else{die "GeneID undefined pattern: ", $_}
	} else{
		warn "no GeneID in line: ", $_;
		sleep 1;
		next;
	}

	# Get ID, replace ID with GeneID:ID
	my $ID = "";
	if (/\tID=([^;]+);/){
		$ID = $1;
		# if (exists $GeneID2ID{$GeneID} ){
		# 	if ($GeneID2ID{$GeneID} eq $ID){next}
		# 	else{die join("\t", $GeneID, $ID, $GeneID2ID{$GeneID}, $_)}
		# }else{
		# 	$GeneID2ID{$GeneID} = $ID;
		# }
		# print $ID;
		# $_ =~ s/\tID=([^;]+);/\tID=$GeneID:$ID;/;

		if (exists $ID2GeneID{$ID} ){
			if ($ID2GeneID{$ID} eq $GeneID){next}
			else{die join("\t", $ID, $GeneID, $ID2GeneID{$ID}, $_)} # check if ID have uniq GeneID, passed
		}else{
			$ID2GeneID{$ID} = $GeneID;
		}
	}else{
		warn "no ID= in line: ", $_;
	}

}

open (OUT, ">ID2GeneID.txt") or die;
foreach my$key (sort keys %ID2GeneID){
	print OUT "$key\t$ID2GeneID{$key}\n";
}
close OUT

