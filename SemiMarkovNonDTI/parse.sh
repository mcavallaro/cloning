
if [ $# -eq 2 ];
then
    prefix=$1
    file=$2

    rm -if $file".txt"


    for i in $(ls $prefix* | sort -nk3 -t"_") ;
    do
        cat $i >> $file"2.txt";
    done

#    python thermodinamicintegration.py $file".txt" $file".pdf"
else
    echo "usage: bash parse.sh [prefix] [out_file_name]"
fi;
