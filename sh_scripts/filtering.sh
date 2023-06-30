#!/bin/bash
## functions
source ./sh_scripts/utils.sh

echo "---- VCF filtering ----"
## read configuration
. parameters.config

echo "Data will be readed from ""$vcf_results_tmp"
echo "and exported in ""$vcffiltpath"
echo "---- ---- ---- ---- ---- ----"

#input path
inputpathlength=`expr length "$vcf_results_tmp"`

#files= `list_files "$vcf_results_tmp" ".vcf"`
####### getting vcf files
vcfsuffix=".vcf"
suffixlength=`expr length "$vcfsuffix"`
untilext=`expr $[suffixlength*-1]`
dirlist=()
for entry in "$vcf_results_tmp"/*
do 
   
   if [ "${entry: $untilext}" == "$vcfsuffix" ];then
      #echo "$entry"
      dirlist+=($entry)
   fi
done

###### go though each file
lengthdirlist=${#dirlist[@]}
for entry in "${dirlist[@]}"
do 
   fnref="${entry%"$vcfsuffix"}"
   newname="$fnref""_flt""$vcfsuffix"
   
   echo "########################################"
   echo "$entry" "$newname"
   echo "########################################"
   vcftools --vcf $entry \
         --recode --remove-indels \
         --maf 0.05 --min-meanDP 10.0 \
         --max-meanDP 50000.0 --minQ 30.0 \
         --min-alleles 2 --max-alleles 2 \
         --max-missing-count 10 -c > $newname
done
