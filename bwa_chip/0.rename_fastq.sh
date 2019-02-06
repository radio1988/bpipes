# first, move raw files to raw

# then sh 0.rename.sh

ls raw/*fastq.gz|sed 's/_.*//'|sort|uniq > sample_name_list

for file in `cat sample_name_list`
do

echo 'For ' $file

ls $file*R1*fastq.gz
echo 'into' ${file/raw\//}_R1.fastq.gz
ln -s $file*R1*fastq.gz  ${file/raw\//}_R1.fastq.gz

ls $file*R2*fastq.gz
echo 'into' ${file/raw\//}_R2.fastq.gz
ln -s $file*R2*fastq.gz  ${file/raw\//}_R2.fastq.gz

echo ' '

done



