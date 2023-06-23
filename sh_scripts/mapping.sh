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
echo "---- Mapping init ----"
## read configuration
. parameters.config

echo "Data will be readed from ""$fastpqualitypath"
echo "and exported in ""$mappingpath"
echo "---- ---- ---- ---- ---- ----"

inputpathlength=`expr length "$fastpqualitypath"`
inputfastextent=`expr length "$fastqextention"`
untilext=`expr $[inputfastextent*-1]`
## output paths check if exists otherwise create
[ ! -d "$mappingpath" ] &&  `mkdir "$mappingpath"`
[ ! -d "$mappingreportpath" ] &&  `mkdir "$mappingreportpath"`
## get genome reference 
genrefpath=`concatenate_paths "$genome_reference_path" "$genome_reference_name"`
echo "$genrefpath"
if [ ! -f "$genrefpath" ];then
    [ ! -f "$genome_reference_name" ] &&  `wget "$genome_url"`
    [ ! -d "$genome_reference_path" ] &&  `mkdir "$genome_reference_path"`
    cp "$genome_reference_name" "$genome_reference_path"
fi
## create genome reference index
[ ! -f "$genrefpath"".ann" ] &&  `bwa-mem2 index "$genrefpath"`

## get files from directory
dirlist=()
for entry in "$fastpqualitypath"/*
do
    if [ "${entry: untilext}" = $fastqextention ];then ## only take those that end with fastq extension
        dirlist+=($entry)
    fi
done
lengthdirlist=${#dirlist[@]}
## find paired files for a specific file and loop the process until there are no more files to process
inputfastextentplus_rq=`expr $[inputfastextent+5]` ## in the cast of fast the names are _R1 so three characters
while [ $lengthdirlist -ne 0 ]
do  
    ## getting paired file names
    pairedfiles=()
    dirlistupt=()
    i=0
    until_r=`expr $[inputfastextentplus_rq*-1]`
    refname="${dirlist[$i]:0:until_r}"
    ## selecting pairs
    for a in "${dirlist[@]}"; do
        if [ "${a:0:until_r}" = "$refname" ];then
            pairedfiles+=($a)
        else
            dirlistupt+=($a)
        fi
    done

    # create names
    samname=()
    samname+=(`concatenate_paths "$mappingpath" "${refname:inputpathlength}.sam"`)
    samname+=(`concatenate_paths "$mappingpath" "${refname:inputpathlength}.bam"`)
    samname+=(`concatenate_paths "$mappingpath" "${refname:inputpathlength}.sorted.bam"`)
    samname+=(`concatenate_paths "$mappingpath" "${refname:inputpathlength}.json"`)
    echo "${samname[0]} ${samname[1]} ${samname[2]} ${samname[3]}"
    echo "$genrefpath ${pairedfiles[0]} ${pairedfiles[1]}"

    # run mapping
    bwa-mem2 mem -t 12 "$genrefpath" "${pairedfiles[0]}" "${pairedfiles[1]}" -M > "${samname[0]}"
    # compress to bam
    samtools view -b "${samname[0]}" > "${samname[1]}"
    # sort and get index
    samtools sort -o "${samname[2]}" "${samname[1]}"
    samtools index "${samname[2]}"
    # report
    samtools flagstat "${samname[2]}" > "${samname[3]}"

    ## remove sam and nor sorted bam files
    [ -e "${samname[0]}" ] && rm "${samname[0]}"
    [ -e "${samname[1]}" ] && rm "${samname[1]}"

    samname="d"

    ## update the paired files
    
    dirlist=()
    for a in "${dirlistupt[@]}"; do
        dirlist+=($a)
    done
    lengthdirlist=${#dirlist[@]}
    echo "################ files left: $lengthdirlist ################"
done 
