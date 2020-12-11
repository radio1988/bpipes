for f in *bedGraph
do 
tail -n +3 $f | cut -f 1,2,3,6 | sort -k1,1 -k2,2n > ${f/bedGraph/bdg} &
done

wait

for f in *bdg
do 
bash bedGraphToBigWig.sh $f ${f/bdg/bw} &
done

wait 

echo "done"
