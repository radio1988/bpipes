cat ../genrich_run3/*narrowPeak| sort -k1,1 -k2,2n|mergeBed>merged.bed
perl bed2gtf.pl merged.bed > merged.gtf
