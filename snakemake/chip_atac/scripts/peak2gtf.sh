for f in *Peak;do perl ../scripts/peak2gtf.pl $f > ${f/_peaks.narrowPeak/.MACS2_BAMPE.gtf};done
