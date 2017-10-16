#!/bin/sh 
# 
# Classifies neighborhoods which include only ony subject gene and multiple
# Copies the into two seperate directories
# 
 
 sources=$1
 singles=$2
 multiples=$3

 for f in $sources*
 do
 	cnt=$(wc -l $f | awk '{print $1}')
 	if [ $cnt -eq 41 ]
 	then
 		cp $f $singles
 	else
 		cp $f $multiples
 	fi
 done