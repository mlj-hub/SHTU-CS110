#!/bin/bash
rm -rf test_record.txt
make fast
make base
for((i=0;i<20;i++))
do 
    echo "----- ${i} run -----" >> test_record.txt
    make fast_test_record
    make base_test_record
    echo "" >> test_record.txt
done