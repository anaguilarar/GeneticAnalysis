#!/bin/bash
## functions
source ./sh_scripts/utils.sh

echo "---- VCF filtering and merging----"
## read configuration
. parameters.config

echo "Data will be readed from ""$vcfmergedpath"
echo "and exported in ""$vcfsummaryresults"
echo "---- ---- ---- ---- ---- ----"

#input path
inputpathlength=`expr length "$vcfmergedpath"`

## output paths
[ ! -d "$vcfsummaryresults" ] &&  `mkdir "$vcfsummaryresults"`

#files= `list_files "$vcf_results_tmp" ".vcf"`
####### getting vcf files
vcfsuffix=".vcf"
suffixlength=`expr length "$vcfsuffix"`
untilext=`expr $[suffixlength*-1]`
dirlist=()
for entry in "$vcfmergedpath"/*
do 
   
   if [ "${entry: $untilext}" == "$vcfsuffix" ];then
      #echo "$entry"
      dirlist+=($entry)
   fi
done

###### go though each file
lengthdirlist=${#dirlist[@]}

sentence=${minmeanDP//,/$'\n'}
for minmeandpval in $sentence
do
   outfiles=()
   for entry in "${dirlist[@]}"
   do 
      fnref="${entry%"$vcfsuffix"}"

      echo "########################################"
      echo "Filtering: ""$entry" 
   
   #vallist= echo $minmeanDP | sed 's/\,//g'

      newname="$fnref""_minDP_""$minmeandpval""_flt""$vcfsuffix"
      new_basename="$(basename $newname)"

      ## output dp paths
      outdpfn=`concatenate_paths "$vcfsummaryresults" "dp_""$minmeandpval"`
      [ ! -d "$outdpfn" ] &&  `mkdir "$outdpfn"`

      outfn=`concatenate_paths "$outdpfn" "$new_basename"`
      echo "Output: ""$outfn"
      echo "########################################"
      vcftools --vcf $entry \
            --recode --remove-indels \
            --maf "$maf" --min-meanDP "$minmeandpval" \
            --max-meanDP 50000.0 --minQ 30.0 \
            --min-alleles 2 --max-alleles 2 \
            --max-missing-count "$mxmissing" -c > $outfn
      outfiles+=($outfn)

   done
   
   fn=""
   for entry in "${outfiles[@]}"
   do
      fn+="-I $entry "
   done

   outfn=`concatenate_paths "$outdpfn" "allchromosomes_minDP_""$minmeandpval""_flt.vcf"`

   echo "########################################"
   echo "Merging: ""$fn" "$outfn"
   java -jar "$picardpath" MergeVcfs \
         $fn -O $outfn
done
