#!/bin/bash
rm -rf fast_test_record.txt
make fast
for((i=0;i<2;i++))
do 
    echo "----- ${i} run -----" >> fast_test_record.txt
    make fast_test_record
    echo "" >> fast_test_record.txt
done