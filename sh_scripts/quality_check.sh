#!/bin/bash
## path
## functions
concatenate_paths() {
   local base_path=${1}
   local sub_path=${2}
   local full_path="${base_path:+$base_path/}$sub_path"
   full_path=$(realpath ${full_path})
   echo $full_path
}
#
echo "---- Quality check init ----"
## read configuration
. parameters.config

echo "Data will be readed from ""$rawfastqpath"
echo "and exported in ""$fastpqualitypath"
echo "---- ---- ---- ---- ---- ----"

inputfastextent=`expr length "$fastqextention"`
untilext=`expr $[inputfastextent*-1]`

#create folders
[ ! -d "$fastpqualitypath" ] &&  `mkdir "$fastpqualitypath"`
#ouputpath_q_r=quality_output_reports
[ ! -d "$fastpqualityreportpath" ] &&  `mkdir "$fastpqualityreportpath"`

## get files from directory
dirlist=()
for entry in "$rawfastqpath"/*
do
    if [ "${entry: untilext}" = $fastqextention ];then ## only take those that end with fastq extension
        dirlist+=($entry)
    fi
done
# check length
inputpathlength=`expr length "$rawfastqpath"`
#pattern length the common pattern for out files is the extension and R1
lengthdirlist=${#dirlist[@]}
## find paired files for a specific file and loop the process until there are no more files to process
inputfastextent=`expr length "$fastqextention"`
inputfastextentplus_r=`expr $[inputfastextent+3]` ## in the cast of fast the names are _R1 so three characters
while [ $lengthdirlist -ne 0 ] ## repeat until there are no more files to process
do  
    pairedfiles=()
    dirlistupt=()
    until_r=`expr $[inputfastextentplus_r*-1]`
    refname="${dirlist[0]:0:until_r}"
    ## selecting pairs
    for a in "${dirlist[@]}"; do
        if [ "${a:0:until_r}" = "$refname" ];then
            pairedfiles+=($a)
        else
            dirlistupt+=($a)
        fi
    done
    ## getting quality names
    ouptuquality=()
    for a in "${pairedfiles[@]}"; do
        #output file names
        untilext=`expr $[inputfastextent*-1]`
        ouptuquality+=(`concatenate_paths "$fastpqualitypath" "${a:inputpathlength: untilext}""_q.fastq.gz"`)
    done
    #output reports
    until_r=`expr $[inputfastextentplus_r*-1]`
    reportnamehtml=`concatenate_paths "$fastpqualityreportpath" "${pairedfiles[0]:inputpathlength: until_r}""_fastp.html"`
    reportnamejson=`concatenate_paths "$fastpqualityreportpath" "${pairedfiles[0]:inputpathlength: until_r}""_fastp.json"`
        
    ## run quality for paired files
    echo "##   Start quality control for: ""${pairedfiles[0]}" "${pairedfiles[1]}"
    fastp -i "${pairedfiles[0]}" -I "${pairedfiles[1]}" \
            -o "${ouptuquality[0]}" -O "${ouptuquality[1]}" \
            -q 25 -n 5 -f 3 -F 3 -h "$reportnamehtml" -j "$reportnamejson"
    
    ## update the paired files
    dirlist=()
    for a in "${dirlistupt[@]}"; do
        dirlist+=($a)
    done
    lengthdirlist=${#dirlist[@]}
    echo "################ files left: $lengthdirlist ################"
    
done 
