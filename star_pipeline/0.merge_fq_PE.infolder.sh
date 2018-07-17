# put biological replicates in folders
# name folders in biologically meaningful ways
# merge R1s and R2s from each folder into Folder_Name_R1/R2.fastq.gz
for dir in */
do
	dir=${dir%/}
	echo $dir 
	echo 'R1s:' $dir/*_R1_*fastq.gz
    cat $dir/*_R1_*fastq.gz > ${dir}_R1.fastq.gz
	echo 'Output: R1' ${dir}_R1.fastq.gz
    
	echo 'R2s:' $dir/*_R2_*fastq.gz
    cat $dir/*_R2_*fastq.gz > ${dir}_R2.fastq.gz
	echo 'Output R2:' ${dir}_R2.fastq.gz
	echo 'done'
	echo ''
done
