#!/usr/bin/env bash

base_dir=/home/hudaiber/Projects/NewSystems/data/cas4/

orgs_file=$base_dir/cas4_orgs.txt

#if [ ! -f $orgs_file ]; then
#    grep -v "#" cas4_census.txt | cut -f7 | sort | uniq > $orgs_file
#fi

####################################################################################
########### Extract annotations for 16s RNA  #######################################
####################################################################################

#for db in Prok1603 GenBank1603 GenMark1603 WGS1603
#do
#    echo
#    echo "Searching 16s RNA in " $db
#    prok_extract_pty -f=4 -p=$db -rty -all $orgs_file | grep '16S' > $base_dir/rty/$db.rty
#    echo "Saved to:"
#    echo $base_dir/rty/$db.rty
#    echo
#done

####################################################################################
########### Retrieve sequences for annotated 16S RNAs  #############################
####################################################################################

#base_dir=/home/hudaiber/Projects/NewSystems/data/cas4/rty/
#
#while read line
#do
#    acc=`echo $line | cut -d" " -f1`
#    coords=`echo $line | cut -d" " -f2 | sed 's/\.\./-/g'`
#
#    _strand=$(echo $line | cut -d" " -f3)
#    if [ "$_strand" == "+" ]; then
#        strand="plus"
#    elif [ "$_strand" == "-" ]; then
#        strand="minus"
#    fi
#
#    org=`echo $line | cut -d" " -f4`
#    contig=`echo $line | cut -d" " -f5`
#
#    echo $org $contig $_strand $coords
#    blastdbcmd -db /panfs/pan1/prokdata/db/all1603.nt -entry $contig -range $coords -strand $strand \
#                > $base_dir/16s_rnas/$org.fna
#
#done < $base_dir/cas4_found_orgs.txt