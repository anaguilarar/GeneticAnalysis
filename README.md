# Variant discovery using GRAS-Di information

The following text aims to explain the process used to identify SNPs and indels in individual samples, resulting in a VCF file containing all the information. The following methodology attempts to replicate the procedure explained by [Takeshima et al. (2022)](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-022-03722-6)

The workflow is as follows:
- Setting up a Conda environment
- Filter low quality reads
- Map to reference
- Call variants per sample
- Consolidate GVCFs across samples
- VCF filtering

## Setting up a Conda environment
The first step was to create an environment containing all the tools used to find variants. The environment was created using conda, an environment management system. All the information needed to install it can be found at: https://docs.conda.io/en/latest/miniconda.html#linux-installers.

Once conda is installed, the next step is to create an environment. Since we will be using the [gatk](https://gatk.broadinstitute.org/hc/en-us) tools to get the VCF files, this requires a set of packages that can be easily installed at the moment of creating a new conda environment. So with that in mind, we will first download the gatk-4.4.0.0 tool [here](https://github.com/broadinstitute/gatk/releases) and then create an environment called

```bash
conda env create -n gatk -f "GATK PATH/"gatkcondaenv.yml
```
where you need to change "GATK PATH" to point to the directory containing the unzipped gatk tool.
For more information on how to install gatk using conda, see the following [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4)

After that, you need to enable the environment and install the latest version of java
```bash
conda activate gatk
sudo apt-get update
sudo apt-get install default-jre
apt install openjdk-17-jdk openjdk-17-jre
```

## Filter low quality reads
To check the quality of the reads, we used the tool called [fastqp](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234), which not only performs quality control, but also allows to apply operations such as filtering, adapter trimming and others.
Installation was done by downloading the package from the conda repository.
```bash
conda install -c bioconda fastp
```
As our purpose was to replicate the above study, the parameters we used were -q 25, -n 5, -f 3, -F 3 and -l 30.
We also exported the quality reports as html and json files. More information can be found in the [fastp documentation](https://github.com/OpenGene/fastp).


```bash
fastp -i paired_file_R1.fastq -I paired_file_R1.fastq \
    -o paired_file_R1_q.fastq -O paired_file_R1_q.fastq \
    -q 25 -n 5 -f 3 -F 3 -h report.html -j report.json
    
```

## Map to reference
Sequence mapping is a process used to detect variantss, such as single-nucleotide polymorphisms (SNPs), indels, etc., between DNA samples and a reference genome. This process was carried out using the tool [BWA-MEM2](https://arxiv.org/pdf/1907.12931.pdf) tool.

To use it, we need to download the package that is available in conda:
```bash
conda install -c bioconda bwa-mem2
```

To perform the sequencing mapping, we need to have the reference genome, so since our study is based on Chinese cabbage, we download the sequence from the Brassicaceae Database (BRAD) project. Specifically, the sequence for Brassica Chiifu V3.5. The following command can be used to obtain it.

```bash
wget http://39.100.233.196:82/download_genome/Brassica_Genome_data/Brara_Chiifu_V3.5/Brara_Chiifu_V3.5.fa.gz
```

Once the reference genome was obtained, the next step was to index the reference using the bwa-mem2 index function.
```bash
bwa-mem2 index Brara_Chiifu_V3.5.fa.gz
```

With the reference genome, the next step is to create a SAM file representing the aligned sequences. For this we will use bwa-mem2 again, but with the mem function, which works using a seeding alignments with maximum exact matches (MEMs) algorithm.

```bash
bwa-mem2 mem -t 12 reference_genome.fa.gz paired_file_R1_q.fastq.gz paired_file_R2_q.fastq.gz > output.sam
```
Where -t is the number of threads the tool will use. For more options please see this [link](https://bio-bwa.sourceforge.net/bwa.shtml).

Since SAM files are large, it is recommended to compress them into a BAM file. For this we will use [samtools](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/pdf/btp352.pdf). This is a tool for processing high throughput sequencing data.

```bash
samtools view -b output.sam > output.bam
```
The next step is to sort and index the BAM files. Sorting by coordinates is used to streamline data processing and avoid loading extra alignments into memory, while indexing function indexes the sorted position in the BAM file [(Li et al. 2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/pdf/btp352.pdf)

```bash
#-o denotes the output filename
samtools sort -o output.sorted.bam output.bam
samtools index output.sorted.bam
```
Finally, a JSON report indicating the mapping percentage is obtained using the flagstat function. available in samtools.
```bash
samtools flagstat output.sorted.bam > output.json
```

## Call variants per sample (GVCFs files per sample)
The identification of variants per sample is done by gatk. To call the variants per sample we use the function HaplotypeCaller. But first it was necessary to add read groups to each bam file using [picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-). The latest version of picard (3.0.0) can be downloaded from this [link](https://broadinstitute.github.io/picard/).

To use the AddOrReplaceReadGroups function, we used the default options mentioned in [gatk](https://gatk.broadinstitute.org/hc/en-us/articles/360035532352-Errors-about-read-group-RG-information).

```bash
java -jar "picard-3.0.0path/picard.jar" AddOrReplaceReadGroups \
            I=bamfile.bam \
            O=bamfile_withreads.bam \
            SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina \
            RGSM=Sample1 \ # this is the sample name, for multiple samples change the name to identify each one
            CREATE_INDEX=True RGPU=unit1

```
The next step was to create a sequence dictionary for the reference genome. This is done by 1) specifying the directory containing the gatk tool 2) decompressing the gs file 3) indexing our reference genome using the CreateSequenceDictionary function 4) creating index.

```bash
#Decompress the gs file
gzip -d Brara_Chiifu_V3.5.fa.gz
gatkpath/gatk CreateSequenceDictionary -R Brara_Chiifu_V3.5.fa
samtools faidx Brara_Chiifu_V3.5.fa
```

Finally, the gvcf file is obtained using the HaplotypeCaller function. We set Java memory allocation to 16 Gb (-Xmx16g). 

```bash
gatkpath/gatk --java-options "-Xmx16g" HaplotypeCaller \
                -R Brara_Chiifu_V3.5.fa \
                -I bamfile_withreads.bam \
                -O output.g.vcf.gs \
                -ERC GVCF \ ## Specify file type
                --native-pair-hmm-threads 10 ## the number of threads 
```

## Consolidate GVCFs across samples

The consolidation of the GVCFS samples was done using the GenomicsDBImport and GenotypeGVCFS functions, both from gatk. The former imports data from one or more samples over at least one genomic interval (chromosome in our case) and results in a GenomicDB datastore, while the latter reads this output and generates a final VCF file. For more information please use this [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035889971--How-to-Consolidate-GVCFs-for-joint-calling-with-GenotypeGVCFs)


```bash
gatkpath/gatk --java-options "-Xmx16g -Xms4g" GenomicsDBImport \
   -V sample1.g.vcf.gs \
   -V sample2.g.vcf.gs \
   -V sample3.g.vcf.gs \
   --genomicsdb-workspace-path vcf_output_path \
   -L A01 ## chromosome reference
```

```bash
gatkpath/gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R Brara_Chiifu_V3.5.fa  \ 
    -V "gendb://vcf_output_path" \
    -O outputfile_A01.vcf 
```

## Filtering
Once the VCF files were obtained per chromosome, the next step was to remove low quality indels, markers with more than N missing SNPs, SNPs with more than two alleles. To accopmlish this, we employed [Vcftools](https://vcftools.sourceforge.net/)

```bash
vcftools --vcf outputfile_A01.vcf \
         --recode --remove-indels \
         --maf 0.05 --min-meanDP 10.0 \
         --max-meanDP 50000.0 --minQ 30.0 \
         --min-alleles 2 --max-alleles 2 \
         --max-missing-count 10 -c > outputfile_A01_flt.vcf
```

## Merging into a single VCF file

```bash
java -jar "picard-3.0.0path/picard.jar" MergeVcfs \
        -I vcftmp/outputfile_A01.vcf \
        -I vcftmp/outputfile_A02.vcf \
        -O output_allchromosomes.vcf
```


## Construction of linkage maps

So far all the steps tat has been appied has been using bash. But the next steps were done in [R](https://www.r-project.org/). 

# How to use bash scritps

The first thing to do is to change the configuration parameters file "parameters.config". These parameters indicate which directory path contains the fastq files, as well as the output directory names. Each parameter is explained below.

`Reference genome`
- *genome_reference_name*: genome reference file name
- *genome_reference_path*: path in which the genome reference is located
- *genome_url*: if the genome reference is not on your local computer, which link can be used to download the data. If you have the data, leave it as "".

`Quality parameters`
- *rawfastqpath*: input path containing the raw fastq files
- *fastpqualitypath*: folder to use for exporting quality files.
- *fastpqualityreportpath*: Folder to use for exporting quality report files.
- *fastqextension*: extension for fastq files. E.g. ".fq.gz".
- ncharactersforparing: Length of the pairing string, e.g. "_R1" is three characters.

`Mapping parameters`
- *mappingpath*: folder used to export bam and bam.bai files
- *mappingreportpath*: folder used to export mapping report files

`GVCF parameters`
- *gatkpath*: path specifying the location of the gatk tool
- *picardpath*: path specifying the location of the picard tool
- *vcf_path*: folder used to export g.vcf.gs files

`VCF merge parameters`
- *vcfmergedpath*: folder to use for exporting vcf files
- *chrn*: total number of chromosomes
- *chrsuffix*: name used to identify the chromosome

Once the parameters has been changed, you can run each bash file, usin bash "bash script.sh". For example for aulity check


```bash
bash sh_scripts/quality_check.sh
```

### Outputs
Inside of each output dir there will be the follogin files
```
├── quality_results
│  ├── ril_1_1_q.fastq
│  ├── ril_1_1_q.fastq
├── quality_reports
│  ├── ril_1.html
│  ├── ril_1.json
```