2019/05/30 summary

## download from NCBI
- dir: /project/umw_nathan_lawson/Rui/zb_genome/ucsc_GRCz11/full/gff_conversion/
- file: GCF_000002035.6_GRCz11_genomic.gff.gz
- Problems:
    - Ugly chr names
    - messy scaffolds, causing problem with Mapping Quality score
    - Onur discovered problems: one transcripts on many chrs (NCBI's problem)

## converted chr_names to UCSC version
- dir: /project/umw_nathan_lawson/Rui/zb_genome/ucsc_GRCz11/full/gff_conversion/
- file: GCF_000002035.6_GRCz11_genomic.ucsc.gff.gz
- code: perl gtf_chrName_change.pl GCF_000002035.6_GRCz11_genomic.gff ucscToRefSeq.txt
- Problems:
    - messy scaffolds

## keep primary chr only
- dir:/project/umw_nathan_lawson/Rui/zb_genome/ucsc_GRCz11/primary/danRer11.Refseq/raw 
- file: GCF_000002035.6_GRCz11_genomic.ucsc.primary.gff.gz
- Problems:
    - GFF not compatible with some downstream code

## convert to GTF with gffread
- dir: /project/umw_nathan_lawson/Rui/zb_genome/ucsc_GRCz11/primary/danRer11.Refseq/raw
- file: CF_000002035.6_GRCz11_genomic.ucsc.primary.gtf.gz
- code: gffread GCF_000002035.6_GRCz11_genomic.ucsc.primary.gff -T -F -O -o GCF_000002035.6_GRCz11_genomic.ucsc.primary.gtf
- can be used for featureCount for gene level quantification, but has some minor issues:
    - transcript_id duplicated, now have both 'rna0' and 'NM_173235.3' in the same line
    - gene_biotype discarded
    - gene_ids become 'gene0' 'gene1'

## fix gene_id
- dir: /project/umw_nathan_lawson/Rui/zb_genome/ucsc_GRCz11/primary/danRer11.Refseq/raw/
- code: gff_GeneID2Name.pl > ID2GeneID.txt
- code: <gtf_fix_gene_id.pl> <ncbi.refseq.gff3.gffread.gtf> <ID2GeneID.txt>

## fix transcript_id
- dir: /nl/umw_nathan_lawson/pub/UCSC_tracks/GRCz11/Improved3pUTRannotations/StartingAnnotations/intermediate/
- code:  zcat GRCz11refSeqUcsc.gtf.gz|perl fix_transcript_id.pl > GRCz11refSeqUcsc.fixed.gtf

## Final version:
- GFF: /project/umw_nathan_lawson/Rui/zb_genome/ucsc_GRCz11/primary/danRer11.Refseq/raw/GCF_000002035.6_GRCz11_genomic.ucsc.primary.gff.gz
    - Problems:
        - should be good, only Onur's problem by NCBI
        - not tested yet
- GTF: /project/umw_nathan_lawson/Rui/zb_genome/ucsc_GRCz11/primary/danRer11.Refseq/raw/
    - Problems:
        - Onur's problem
        - no biotype


