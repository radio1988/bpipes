todo: 
- merge peaks then DESeq2
- t_vs_c step has to take 1 rep only
- MEME summit 

Notes:
- sample names [a-Z0-9_.] only
- peaks: called by macs2, narrow, broad
- clean_peaks: filtered with blacklist
- clean.real.peaks: filtered with CPM, Pulldown > IgG
- final_anno.xlsx: clean.real.peaks annotated by GTF
	- distanceToSite: distanceToTSS
	- Norm: DESeq2 normed count, size_factor based on total reads in peaks
