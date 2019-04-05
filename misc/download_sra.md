## Experiences
- If you download with fastq-dump, do not initiate more than 3 jobs at a time (2019), or some error likely to show up.
- `fastq-dump` takes large space in $HOME/ncbi, you need to create softlink to borrow space from a larger storage
- Using ENA is preferred, compared with fastq-dump: 
  - no sra.cache files
  - more robust?
- Trick: with 2+ runs for a sample, the scripts just don't work, e.g. SRS454725
  

