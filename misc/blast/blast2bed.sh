# gffread -w gencode.v29.primary_assembly.annotation.fixed.fa -g hg38.primary.fa gencode.v29.primary_assembly.annotation.fixed.gtf

# makeblastdb -in Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -input_type fasta -dbtype nucl -title GRCh38

#blastn -task blastn -word_size 7 -num_threads 4 -query primers.fa -db /Users/rui/genome/human/hg38_primary_ucsc/hg38.primary.fa -outfmt 11  > primers.blast.asn
blastn -task blastn -word_size 7 -num_threads 4 -query primers.fa -db /Users/rui/genome/human/hg38_primary_ucsc/hg38.primary.fa -outfmt 7  > primers.blast.txt
blastn -task blastn -word_size 7 -num_threads 4 -query primers.fa -db /Users/rui/genome/human/hg38_primary_ucsc/hg38.primary.fa -outfmt 17 > primers.blast.sam

perl blast7_2_bed.pl primers.blast.txt > primers.blast.bed
sort -k1,1 -k2,2n primers.blast.bed > temp.bed && mv temp.bed primers.blast.bed

bedtools intersect -a primers.blast.bed -b ~/genome/human/hg38_primary_ucsc/gencode.v29.primary_assembly.annotation.fixed.exon.gtf -wa > primers.blast.exon.bed 
