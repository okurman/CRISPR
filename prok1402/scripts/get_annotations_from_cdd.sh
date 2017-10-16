#!/bin/sh 
#
# Get annotation lines from Prok1402.*.ccp.lst files which 
# correspond to profile identifiers read from the file given as
# argument. Which is obtained with:
# cat cdd_* | awk '{print $2}' | sort | uniq > profile_ids.txt
#

if [ $# -eq 0 ]; then
	echo "Provide paths to 1: Archive of annotations 2: list of profile ids to search annotations for."
	exit 1
fi

archive_file=$1
profile_ids=$2

if [[ -d ./tmp ]]; then
	rm -rf tmp
fi

mkdir tmp

while read -r line
do
	cur_id=$line
	tmp_file=$cur_id".tmp"
	grep $cur_id $archive_file > ./tmp/$tmp_file
	echo $tmp_file
done < $profile_ids

cat tmp/*.tmp > selected_annotations.csv
rm -rf tmp