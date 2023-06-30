#!/bin/bash
## functions
concatenate_paths() {
   local base_path=${1}
   local sub_path=${2}
   local full_path="${base_path:+$base_path/}$sub_path"
   full_path=$(realpath ${full_path})
   echo $full_path
}

list_files(){
   local base_path=${1}
   local suffix=${2}

   suffixlength=`expr length "$suffix"`
   dirlist=()
   untilext=`expr $[suffixlength*-1]`
   #echo "$untilext"
   #echo "$base_path"
   files="${1}"/*
   lengthdirlist=${#files[@]}
   echo $lengthdirlist
   for entry in "${files[@]}"
   do 
      #echo "$entry"
      if [ "${entry: $untilext}" == "$suffix" ];then
         
         dirlist+=($entry)
      fi
   done
   echo $dirlist
}