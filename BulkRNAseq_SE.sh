#!/bin/bash

# Exit on error
set -e

## Download
# Set variable listing all desired accessions
accessions=(SRR497699 SRR49770{0..9})

# Download raw fastqs
for accession in "${accessions[@]}"; do
echo "Processing accession: $accession"
fastq-dump --gzip $accession
done

## QC
# FastQC
mkdir ./fastqc_initial
export PATH=/programs/FastQC-0.12.1:$PATH

fastqc -t 2 *.fastq* -o ./fastqc_initial

# Compile into a single report with MultiQC
export PYTHONPATH=/programs/multiqc-1.15/lib64/python3.9/site-packages:/programs/multiqc-1.15/lib/python3.9/site-packages
export PATH=/programs/multiqc-1.15/bin:$PATH

multiqc .

## Trim reads with Trimmomatic
# Set variables
input_dir="/workdir/sml327/PATH"
# Get unique SRR-style sample IDs from filenames in the directory
donors=($(ls "$input_dir" | grep -oP 'SRR\d+' | sort -u))

# Export Trimmomatic paths for use inside GNU parallel
export TRIM_JAR=/programs/trimmomatic/trimmomatic-0.39.jar
export TRIM_ADAPTERS=/programs/trimmomatic/adapters/TruSeq3-SE.fa

# Run in parallel (15 jobs)
printf "%s\n" "${donors[@]}" | parallel -j 15 '
  sample={}
  echo "Processing sample: ${sample}"

  in_fastq="'"$input_dir"'/${sample}.fastq"
  out_fastq="'"$input_dir"'/${sample}_trimmed.fastq"

  echo "Input:  ${in_fastq}"
  echo "Output: ${out_fastq}"

  java -jar "${TRIM_JAR}" SE -phred33 \
    "${in_fastq}" "${out_fastq}" \
    ILLUMINACLIP:"${TRIM_ADAPTERS}":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
'

## Set up SQuIRE
cat $(dirname $CONDA_PREFIX)/.condarc
#Output should be:
#channels:
# - conda-forge
# - bioconda
# - defaults
conda activate squire

# Check versions to make sure they're compatible with SQuIRE
R --version && python --version && stringtie --version && bedtools --version && samtools --version && STAR --version
#Output should be:
#R 3.4.1 / Python 2.7 / Stringtie 1.3.3 / Bedtools 2.25.0 / Samtools 1.1 / STAR 2.5.3a

## SQuIRE Map to map the reads to hg38
# Fetch/clean directories are from a previous SQuIRE run
# Set variables
annotation="/workdir/sml327/SQuIRE/squire_fetch/hg38_refGene.gtf"
fetch_dir="/workdir/sml327/SQuIRE/squire_fetch"
# Variables below will need adjustment
input_dir="/workdir/sml327/PATH"
output_dir="/workdir/sml327/PATH/squire"
read_length="90"
threads_per_job="4"
num_parallel="2"

# Run in parallel 
parallel -j $num_parallel '
  echo "Processing file: {}"
  squire Map -1 {} \
    -o '"$output_dir"' \
    -f '"$fetch_dir"' \
    -r '"$read_length"' \
    -b hg38 \
    -g '"$annotation"' \
    -p '"$threads_per_job"' -v
' ::: "$input_dir"/*_trimmed.fastq

mkdir ./trim
mv *trimmed* ./trim
cd ./trim

## Determine strandedness
export PYTHONPATH=/programs/RSeQC-5.0.4/lib64/python3.9/site-packages:/programs/RSeQC-5.0.4/lib/python3.9/site-packages
export PATH=/programs/RSeQC-5.0.4/bin:$PATH
file="/workdir/sml327/PATH/squire/SAMPLE_trimmed.bam"

infer_experiment.py -r /workdir/sml327/SQuIRE/squire_fetch/hg38_refGene.bed -i "$file"

## SQuIRE Count to count gene and TE reads
# Set variables
clean_dir="/workdir/sml327/SQuIRE/squire_clean"
fetch_dir="/workdir/sml327/SQuIRE/squire_fetch"
# Variables below will need adjustment
map_dir="/workdir/sml327/PATH/squire/squire_map"
read_length="90"
threads="80"
strand="2"

for file in "$map_dir"/*_trimmed.bam; do
  # Extract the base name without extension or path
  base=$(basename "$file" "_trimmed.bam")
  echo "Counting for sample: ${base}"

  squire Count \
    -m "$map_dir" \
    -c "$clean_dir" \
    -f "$fetch_dir" \
    -r "$read_length" \
    -n "${base}" \
    -b hg38 \
    -p "$threads" \
    -s "$strand" \
    -v
done

## SQuIRE Call for DEA
# Set variables
count_dir="/workdir/sml327/PATH/squire/squire_map/squire_count"
group1="basename1,basename2,basename3"
group2="basename4,basename5,basename6"
condition1="C1"
condition2="C2"
projname="NAME"

# Locus level
squire Call \
  -1 "$group1" \
  -2 "$group2" \
  -A "$condition1" \
  -B "$condition2" \
  -i "$count_dir" \
  -p 40 \
  -N "$projname" \
  -f pdf \
  -t TRUE \
  -v

# Subfamily level
squire Call \
  -1 "$group1" \
  -2 "$group2" \
  -A "$condition1" \
  -B "$condition2" \
  -i "$count_dir" \
  -s True \
  -p 40 \
  -N "$projname" \
  -f pdf \
  -v

echo "RNA-seq pipeline completed."
