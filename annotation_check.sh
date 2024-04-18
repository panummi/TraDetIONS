#!/bin/bash

bed=$1
vcf=$2
txt=$3
tefile=$4

while IFS= read -r id
	do
	origte=$(grep -w $id <(less $vcf ) | cut -f3,8 | sed 's/SVclass=/\t/g' | cut -f1,3 | cut -f1 -d';')
	echo $origte | tr '\n' '\t'
	te=$(echo $origte | awk '{print $2}')
    firstmatch=$(grep -w "$te" $tefile | cut -f2)
	[[ ! -z "$firstmatch" ]] && echo "$firstmatch" | tr '\n' '\t' || echo "-" | tr '\n' '\t'
	transp=$(grep -w $id $bed | grep transposon | cut -f7- | cut -f1 -d':')
	[[ ! -z "$transp" ]] && echo "$transp" | tr '\n' '\t' || echo "-" | tr '\n' '\t'
	matchte=$(grep -w "$transp" $tefile | cut -f2)
	[[ ! -z "$matchte" ]] && echo "$matchte" | tr '\n' '\t' || echo "-" | tr '\n' '\t'
	if [[ $firstmatch == $matchte ]];then
		echo "Match" | tr '\n' '\t'
	else
		echo "No match" | tr '\n' '\t'
	fi
	poly=$(grep -w $id $bed | grep poly -c)
	en=$(grep -w $id $bed | grep endonuclease -c)
	tsds=$(grep -w $id $bed | grep target_site | cut -f7 | cut -f2 -d':' | sed 's/-/\t/g' | awk '{print $2-$1}')
	array=($tsds)
	propertsd=0
	for element in "${array[@]}"
	do
    		if [[ $element -gt 4 && $element -lt 26 ]];then
			((propertsd=propertsd+1))
		fi
	done
	if [[ $en -gt 0 || $poly -gt 0 || $propertsd -gt 0 ]];then
		echo "hallmarks" | tr '\n' '\t'
	else
		echo "no_hallmarks" | tr '\n' '\t'
	fi
	echo
done < $txt
