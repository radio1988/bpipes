for f in ../*gz;do echo $f; zcat $f|head -n 4000000 |gzip -c > ${f/\./};done
