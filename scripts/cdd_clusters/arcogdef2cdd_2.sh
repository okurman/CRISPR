#!/bin/sh

arcog_def_file=~/data/Archea/arCOG/ar14.arCOGdef.tab
funclass_file=~/data/Archea/arCOG/ar14/funclass.tab


for class_code in $(cut -f1 $funclass_file);
do
	# class_code=$1
	# echo $class_code
	out_fname=$class_code'_arcog_desc.txt';

	count=`grep -cE '^\w+\t'$class_code $arcog_def_file`
	if [ $count -ne 0 ]
	then
		echo $class_code, $count
		grep -E '^\w+\t'$class_code $arcog_def_file | \
			awk -F"\t" '{print $4}' | \
			sed 's/[U,u]ncharacterized//g' | \
			sort | uniq > tmp/$out_fname
	fi

done;

#cat tmp/*_arcog_desc.txt > cdd_to_class.tab;
#rm tmp/*_arcog_desc.txt;
