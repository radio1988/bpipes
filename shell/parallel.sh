i=0
for f in *bdg
do
let "i++"

echo sort -k1,1 -k2,2n $f > ${f/bdg/sort.bdg}
sort -k1,1 -k2,2n $f > ${f/bdg/sort.bdg} &

if !(($i % 4)); then
wait
fi

done
