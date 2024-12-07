

if [ $# -eq 1 ]
then
rm  $1"aggr.txt"
for f in $(ls  $1* | sort -nk 5 -t _)
do
   tail $f -n 1 >> $1"aggr.txt"
done

fi
