i=3
orderl='0,1.2'
orderr=',2.2'
for k in $(ls *.txt)
do
    if [ -a final.results ]
    then
            join -a1 -a2 -e "0" -o "$orderl$orderr" final.results $k  > tmp.res
            orderl="$orderl,1.$i"
            i=$((i+1))
            mv tmp.res final.results
    else
            cp $k final.results
    fi
done
