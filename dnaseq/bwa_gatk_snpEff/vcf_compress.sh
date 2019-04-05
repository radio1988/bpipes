mkdir -p compressed 
for file in *vcf; do bgzip -c $file > ./compressed/${file}.gz;done
for file in compressed/*gz; do tabix -p vcf $file;done
