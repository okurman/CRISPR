#!/bin/sh 
# 
# Get neighborhoods from pty archive file, for given GIs 
# 
 
if [ $# -eq 0 ]; then
	echo "Usage: get_neighborhoods.sh annot_file flank_length out_path"
	echo "Input the list of annotations for which to retrieve the neighborhoods, size of flanking regions"
	echo "and a path to save the neighborhood files for each given GI"
	echo "As source of pty, uses respective directory from  /Users/hudaiber/data/Pty/genomes"
	echo "Format of annotation file must be of standard COG annotation"
	exit 1
fi

pty_dir='/Users/hudaiber/data/Pty/genomes/'

annot_file=$1
flank=$2
out_path=$3

touch found.log
touch not_found.log

cnt=1

cut -d , -f1,2 $annot_file | sort | uniq | while read line
do
	cnt=$((cnt+1))
	
	# for debugging
	# if [[ $cnt -lt 9237 ]]; then
	# 	continue
	# fi

	gi=$(echo $line | cut -d, -f1)
	org_name=$(echo $line | cut -d, -f2)
	# echo $org_name
	# echo $gi
	# echo $pty_dir$org_name/
	# exit 0
	if [ -d "$pty_dir$org_name/" ]; then
		neigh_file=$out_path$gi'.pty'
		scr_name=$(grep -l $gi $pty_dir$org_name/*.pty)
		grep -B $flank -A $flank $gi $scr_name > $neigh_file
		echo $cnt $gi $org_name >> found.log
		# exit 0
	else
		echo $cnt $gi $org_name >> not_found.log
	fi
done