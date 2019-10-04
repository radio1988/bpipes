echo "usage: sort_gtf.sh  name.gtf"
echo "output into:" ${1/gtf/s.gtf} "; then, rename back to $1"
sleep 5s
cat $1 | sort -k1,1 -k4,4n > ${1/gtf/s.gtf} && mv ${1/gtf/s.gtf} $1
