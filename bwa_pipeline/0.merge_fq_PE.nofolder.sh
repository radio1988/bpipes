ls *fastq.gz|sed 's/_.*//'|sort|uniq > sample_name_list

for file in `cat sample_name_list`
do
	echo 'For ' $file

	ls $file*R1*fastq.gz
	cat $file*R1*fastq.gz > ${file}_R1.fastq.gz
	echo 'into' ${file}_R1.fastq.gz

	ls $file*R2*fastq.gz
	cat $file*R2*fastq.gz > ${file}_R2.fastq.gz
	echo 'into' ${file}_R2.fastq.gz

	echo ' '
done


mkdir -p raw
mv *L00*_00* raw
