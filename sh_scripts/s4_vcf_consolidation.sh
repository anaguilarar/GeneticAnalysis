#!/bin/bash
## functions
source ./sh_scripts/utils.sh

echo "---- GVCF consolidation ----"
## read configuration
. parameters.config

echo "Data will be readed from ""$vcf_path"
echo "and exported in ""$vcfmergedpath"
echo "---- ---- ---- ---- ---- ----"

#input path
inputpathlength=`expr length "$vcf_path"`
### setting reference fasta name name
genrefpath=`concatenate_paths "$genome_reference_path" "$genome_reference_name"`
genreference_basename="$(basename $genrefpath)"
gen_reference_fasta=`concatenate_paths "$genome_reference_path" "fasta"`
newrefpath="$(concatenate_paths "$gen_reference_fasta" "$genreference_basename")"
newrefpathfa="${newrefpath%".gz"}"

#output path
[ ! -d "$vcfmergedpath" ] &&  `mkdir "$vcfmergedpath"`

### list files
dirlist=()
for entry in "$vcf_path"/*
do
   if [ "${entry: -9}" == ".g.vcf.gz" ];then
      dirlist+=($entry)
   fi
done

#### creating samples merge command
fn=""
samplename=()
for a in "${dirlist[@]}"; do
   fn+=" -V $a "
   ## correcting sample names
   #creating new name for sample
   sringli=`expr $[inputpathlength+1]`
   samplename="S_${a:$sringli: -9}"
   vcfname="${a%".gz"}"
   gzip -d -k $a
   #checking if the reference sample name is in the file: if it has the same name it will be chjanged by a new one
   val=`grep -R "$samplenameref" $vcfname`
   val=`expr length "$val"`
   if [ $val -ge 1 ];then
      echo "Correcting sample names for: ""$a"
      sed -i "s/$samplenameref/$samplename/g" $vcfname
      `rm "$a"`
      bgzip $vcfname $a
      #removing previous index
      indexname="$a"".tbi"
      [ -f "$indexname" ] &&  `rm "$indexname"`
      gatk-4.4.0.0/gatk-4.4.0.0/gatk IndexFeatureFile -I $a
   else
      `rm "$vcfname"`
   fi

done

#

# chromosomes

for (( i=1; i<=$chrn; i++ )); do
   if [ $i -le 9 ];then
      crh="A0$i"
   else
      crh="A$i"
   fi
   
   vcf_output="vcf_$crh"
   [ -d "$vcf_output" ] &&  `rm "$vcf_output"` #outputpath must be a new folder

   outfile=`concatenate_paths "$vcfmergedpath" "rils_$crh.vcf"`

   echo "########################################"
   echo "$crh" "$vcf_output" "$outfile"
   echo "########################################"
   
   ##### gatk
   gatk-4.4.0.0/gatk-4.4.0.0/gatk --java-options "-Xmx16g -Xms4g" GenomicsDBImport \
      $fn  --genomicsdb-workspace-path $vcf_output -L $crh
   
   echo "the vcf will be saved as: $outfile"
   gatk-4.4.0.0/gatk-4.4.0.0/gatk --java-options "-Xmx16g" GenotypeGVCFs -R $newrefpathfa \
      -V "gendb://$vcf_output" -O "$outfile"

done




