for f in *deseq.csv; do perl filter_sig.pl $f > ${f/csv/sig.csv};done
