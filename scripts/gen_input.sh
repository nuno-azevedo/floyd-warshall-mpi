#!/bin/sh

echo $1 > Input$1.txt
for i in $(seq 1 $1); do
    shuf -r -i 0-50 -n $1 | xargs >> Input$1.txt
done
