singularity exec ~/singularity/hand_sandbox faToTwoBit hg38.primary.cap.fa hg38.primary.cap.2bit
cat hg38.primary.fa|perl -lane 'if (/^>/){print;next}; print uc $_' > hg38.primary.cap.fa
# https://genome.ucsc.edu/goldenpath/help/twoBit.html
