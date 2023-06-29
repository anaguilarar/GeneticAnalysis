#!/bin/bash
## functions
concatenate_paths() {
   local base_path=${1}
   local sub_path=${2}
   local full_path="${base_path:+$base_path/}$sub_path"
   full_path=$(realpath ${full_path})
   echo $full_path
}
#
echo "---- GVCF creation per sample init ----"
## read configuration
. parameters.config

echo "Data will be readed from ""$mappingpath"
echo "and exported in ""$vcf_path"
echo "---- ---- ---- ---- ---- ---- --- ----"


inputpathlength=`expr length "$mappingpath"`
# genome referece 
genrefpath=`concatenate_paths "$genome_reference_path" "$genome_reference_name"`
genreference_basename="$(basename $genrefpath)"
## output paths
[ ! -d "$vcf_path" ] &&  `mkdir "$vcf_path"`
## gen fasta
gen_reference_fasta=`concatenate_paths "$genome_reference_path" "fasta"`
[ ! -d "$gen_reference_fasta" ] &&  `mkdir "$gen_reference_fasta"`

## picard
######## SETTING REFERENCE
# copy file
sringli=`expr length "$genreference_basename"`
newrefpath="$(concatenate_paths "$gen_reference_fasta" "$genreference_basename")"
newrefpathfa="${newrefpath%".gz"}"
# copy and decompressing the reference genome
[[ ! -f "$newrefpathfa" ]] &&  `cp "$genrefpath" "$gen_reference_fasta"` &&  `gzip -d "$newrefpath"`
echo "${newrefpathfa%".fa"}"
# creating dict
[[ ! -f "${newrefpathfa%".fa"}"".dic" ]] &&  `$gatkpath CreateSequenceDictionary -R "$newrefpathfa"`
# creating index
[[ ! -f "$newrefpathfa"".fai" ]] && `samtools faidx "$newrefpathfa"`

####### getting bam files
dirlist=()
inputbamextent=`expr length "$bamextension"`
untilext=`expr $[inputbamextent*-1]`

for entry in "$mappingpath"/*
do
    if [ "${entry: untilext}" = $bamextension ];then ## only take those that end with fastq extension
        dirlist+=($entry)
    fi
done
###### go though each file
lengthdirlist=${#dirlist[@]}
while [ $lengthdirlist -ne 0 ]
do 
    i=0
    pairedfiles=()
    dirlistupt=()
    for a in "${dirlist[@]}"; do
        fnref="${dirlist[$i]%".sorted.bam"}"
        fni="${a%".sorted.bam"}"
        
        if [ $fni = $fnref ];then
            pairedfiles+=($a)
        else
            dirlistupt+=($a)
        fi
    done

    ## setting names
    # reads
    fnreads="${pairedfiles[0]%".sorted.bam"}"".sorted_rg.bam"

    vcffile="$(basename "${pairedfiles[0]}")"
    fnvcf="${vcffile%".sorted.bam"}"".g.vcf.gz"
    fnvcf="$(concatenate_paths "$vcf_path" "$fnvcf")"
    echo "########################################"
    echo "${pairedfiles[0]}" "$fnreads" "$fnvcf"
    echo "########################################"

    samplename="S_""${vcffile%".sorted.bam"}"
    ##### picard
    sringli=`expr $[inputpathlength+1]`
    samplename="${a:$sringli:$untilext}"
    java -jar "$picardpath" AddOrReplaceReadGroups \
            -I "${pairedfiles[0]}" \
            -O "$fnreads" \
            -SORT_ORDER coordinate -RGID foo -RGLB bar -RGPL illumina \
            -RGSM "$samplename" -CREATE_INDEX True -RGPU unit1

    ##### gatk
    ulimit -c unlimited
    $gatkpath --java-options "-Xmx16g" HaplotypeCaller \
                -R "$newrefpathfa" \
                -I "$fnreads" \
                -O "$fnvcf" \
                -ERC GVCF \
                --native-pair-hmm-threads 10

    [ -e "$fnreads" ] && rm "$fnreads"
    [ -e "${fnreads%".bam"}"".bai" ] && rm "${fnreads%".bam"}"".bai"

    dirlist=()
    for a in "${dirlistupt[@]}"; do
        fni2="${a%".sorted.bam.bai"}"
        if [[ $fni2 != $fnref ]];then
            dirlist+=($a)
        fi
    done
    lengthdirlist=${#dirlist[@]}
    echo "################ left: $lengthdirlist ################"
done 