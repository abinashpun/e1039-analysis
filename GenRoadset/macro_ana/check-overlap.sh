#!/bin/bash

RS1=${1:-'123'}
RS2=${2:-'129'}
echo "$RS1 vs $RS2"

for GRP in pos_top pos_bot neg_top neg_bot ; do
    FN1=roadset/$RS1/$GRP.txt
    FN2=roadset/$RS2/$GRP.txt
    N1=$(cat $FN1 | wc -l)
    N2=$(cat $FN2 | wc -l)
    N3=$(grep -f <(awk '{print $1;}' $FN1 | sort) <(awk '{print $1;}' $FN2 | sort) | wc -l)
    P1=$(( 100 * N3 / N1))
    P2=$(( 100 * N3 / N2))
    echo "$GRP  $N1  $N2  $N3  $P1  $P2"
done
