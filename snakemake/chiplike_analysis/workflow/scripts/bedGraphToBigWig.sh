module load singularity/singularity-current
bedGraphToBigWig="singularity exec $HOME/singularity/hand_sandbox bedGraphToBigWig"
SizeFile=$HOME/mccb/genome/Mus_musculus_UCSC_mm10/mm10.chrom.sizes.txt
input=$1
output=$2

$bedGraphToBigWig $input $SizeFile $output
