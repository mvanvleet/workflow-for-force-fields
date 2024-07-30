#!/bin/bash
for file in $@
do
    #molecule=`dirname $file | tr 'a-z' 'A-Z'`
    homo=`grep HOMO $file | head -1 | awk '{print $3}'`
    ip=`grep -m 1 "SETTING E_IONIZATION" $file | awk '{print $4}'`


    printf "%14.6f %16.8f\n" "$homo" "$ip"
done
