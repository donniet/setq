#!/bin/bash

f=${1:-set_counter_collection.csv}

echo "    a~b,      a,      b,      x,     aa,     ab,     ax,     ba,     bb,     bx,     xa,     xb,     xx,total" > ${f}

for i in $(seq 0 256);
do
    ratio=`awk -v var1=${i} -v var2=256 'BEGIN { print  ( var1 / var2 ) }'`
    ./set_counter abx 43 ${ratio} >> ${f}
done

echo "completed."