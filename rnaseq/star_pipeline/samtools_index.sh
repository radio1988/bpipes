for file in *bam
do
samtools index $file $file.bai &
done
