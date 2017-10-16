#!/bin/sh

arcog_def_file=~/data/Archea/arCOG/ar14.arCOGdef.tab
funclass_file=~/data/Archea/arCOG/ar14/funclass.tab


for class_code in $(cut -f1 $funclass_file);
do
	# echo $class_code
	out_fname=$class_code'_cdd_desc.txt';

	grep -E '^\w+\t'$class_code $arcog_def_file | \
			awk -F"\t" '{print $6" "$7}' |  \
			tr ',' '\n' | tr ' ' '\n' | grep -v '^$' | \
			sed 's/$/	'$class_code'/' > $out_fname;

	grep -E '^\w+\t'$class_code $arcog_def_file | \
			awk -F"\t" '{print $5}' |  grep '^COG' | \
			awk '{print substr($0,1,3)substr($0,5)}' | \
			sed 's/$/	'$class_code'/' >> $out_fname

done;

cat *_cdd_desc.txt > cdd_to_class.tab;
rm *_cdd_desc.txt;