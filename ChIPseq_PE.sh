#!/bin/bash

# Exit on error
set -e

## Download
# Set variable listing all desired accessions
accessions=(SRR497699 SRR49770{0..9})

# Download raw fastqs
for accession in "${accessions[@]}"; do
  echo "Processing accession: $accession"
  fastq-dump --split-files --gzip "$accession"
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
mv multiqc_report.html GSEID_report.html  # edit GSEID

cd ..

## Trim reads with Trimmomatic

echo "Trimming paired-end reads"

# Simple loop
for file in *_1.fastq; do
  [ -e "$file" ] || continue
  base=$(basename "$file" _1.fastq)
  echo "Processing pair: ${base}_1.fastq and ${base}_2.fastq"

  java -jar /programs/trimmomatic/trimmomatic-0.39.jar PE -phred33 \
    "${base}_1.fastq" "${base}_2.fastq" \
    "${base}_1_paired_trimmed.fastq" "${base}_1_unpaired_trimmed.fastq" \
    "${base}_2_paired_trimmed.fastq" "${base}_2_unpaired_trimmed.fastq" \
    ILLUMINACLIP:/programs/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# OR: GNU parallel version
#: <<'PE_TRIM_PARALLEL'
#ls *_1.fastq | sed 's/_1.fastq//' | parallel -j 4 '
#  echo "Processing {1}_1.fastq and {1}_2.fastq" ;
#  java -jar /programs/trimmomatic/trimmomatic-0.39.jar PE -threads 4 -phred33 \
#    {1}_1.fastq {1}_2.fastq \
#    {1}_1_paired_trimmed.fastq {1}_1_unpaired_trimmed.fastq \
#    {1}_2_paired_trimmed.fastq {1}_2_unpaired_trimmed.fastq \
#    ILLUMINACLIP:/programs/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
#    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#'
#PE_TRIM_PARALLEL

mkdir -p ./trim
mv *trimmed* ./trim
cd ./trim

## Align reads with Bowtie2

echo "Aligning paired-end reads"

index="/workdir/sml327/Homo_Genome/Human_BAM_Index/hg38_bt2/GRCh38_noalt_as"
input_dir="/local/workdir/sml327/PATH/trim"
outdir="${input_dir}/bam"
mkdir -p "$outdir"

export PATH=/programs/bowtie2-2.5.1-linux-x86_64:$PATH
SAMTOOLS=/programs/samtools-1.9-r9/bin/samtools

BOWTIE2_THREADS=40
SORT_THREADS=4

for r1 in "$input_dir"/*_1_paired_trimmed.fastq*; do
  [ -e "$r1" ] || continue
  r2="${r1/_1_paired_trimmed/_2_paired_trimmed}"

  if [ ! -f "$r2" ]; then
    echo "WARNING: Mate pair for $r1 not found, skipping."
    continue
  fi

  base=$(basename "$r1")
  sample="${base%%.*}"
  sample="${sample%%_1_paired_trimmed}"

  out_sorted="${outdir}/${sample}_sorted.bam"
  out_filtered="${outdir}/${sample}_filtered.bam"
  log_file="${outdir}/${sample}_bowtie2.log"

  echo "Processing $sample"

  bowtie2 --phred33 -p "$BOWTIE2_THREADS" --very-sensitive-local \
    -x "$index" -1 "$r1" -2 "$r2" \
    --rg-id "$sample" \
    --rg SM:"$sample" \
    --rg PL:ILLUMINA \
    --rg LB:lib1 \
    2> "$log_file" \
  | $SAMTOOLS view -bS - \
  | $SAMTOOLS sort -@ "$SORT_THREADS" -o "$out_sorted" -

  $SAMTOOLS index "$out_sorted"

  $SAMTOOLS view -h "$out_sorted" \
    | awk '$1 ~ /^@/ || $3 != "chrM"' \
    | $SAMTOOLS view -b -o "$out_filtered" -

  $SAMTOOLS sort -@ "$SORT_THREADS" -o "$out_filtered" "$out_filtered"
  $SAMTOOLS index "$out_filtered"

  echo "Finished $sample: ${out_filtered}"
done

# OR: GNU parallel version
#: <<'PE_BOWTIE_PARALLEL'
#input_dir="/workdir/sml327/PATH/trim"
#outdir="${input_dir}/bam"
#mkdir -p "$outdir"
#
#export PATH=/programs/bowtie2-2.5.1-linux-x86_64:$PATH
#SAMTOOLS=/programs/samtools-1.9-r9/bin/samtools
#
#BOWTIE2_THREADS=4
#SORT_THREADS=2
#
#ls "$input_dir"/*_1_paired_trimmed.fastq* \
#  | sed 's/_1_paired_trimmed.fastq.*//' \
#  | parallel -j 4 --eta '
#    sample=$(basename {})
#    r1={}_1_paired_trimmed.fastq
#    r2={}_2_paired_trimmed.fastq
#
#    out_sorted="'$outdir'/${sample}_sorted.bam"
#    out_filtered="'$outdir'/${sample}_filtered.bam"
#    log_file="'$outdir'/${sample}_bowtie2.log"
#
#    echo "Processing $sample"
#
#    bowtie2 --phred33 -p '$BOWTIE2_THREADS' --very-sensitive-local \
#      -x '$index' -1 "$r1" -2 "$r2" \
#      --rg-id "$sample" \
#      --rg SM:"$sample" \
#      --rg PL:ILLUMINA \
#      --rg LB:lib1 \
#      2> "$log_file" \
#    | '$SAMTOOLS' view -bS - \
#    | '$SAMTOOLS' sort -@ '$SORT_THREADS' -o "$out_sorted" -
#
#    '$SAMTOOLS' index "$out_sorted"
#
#    '$SAMTOOLS' view -h "$out_sorted" \
#      | awk "\$1 ~ /^@/ || \$3 != \"chrM\"" \
#      | '$SAMTOOLS' view -b -o "$out_filtered" -
#
#    '$SAMTOOLS' sort -@ '$SORT_THREADS' -o "$out_filtered" "$out_filtered"
#    '$SAMTOOLS' index "$out_filtered"
#
#    echo "Finished $sample: ${out_filtered}"
#  '
#PE_BOWTIE_PARALLEL

## Deduplicate with Picard MarkDuplicates

echo "Deduplicating"

input_dir="$(pwd)"
outdir="${input_dir}/dedup"
mkdir -p "$outdir"

for bam in ${input_dir}/*.bam; do
  [ -e "$bam" ] || continue
  base=$(basename "$bam" .bam)

  echo "Deduplicating $base"

  java -jar /programs/picard-tools-2.26.1/picard.jar MarkDuplicates \
    I="$bam" \
    O="${outdir}/${base}_dedup.bam" \
    M="${outdir}/${base}_dup_metrics.txt" \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT

  /programs/samtools-1.9-r9/bin/samtools sort -@ 4 -o "${outdir}/${base}_dedup.sorted.bam" "${outdir}/${base}_dedup.bam"
  /programs/samtools-1.9-r9/bin/samtools index "${outdir}/${base}_dedup.sorted.bam"

  echo "Finished $base"
done

## Post-alignment QC

echo "Running samtools QC (flagstat/idxstats)"

input_dir="/local/workdir/sml327/PATH/dedup"
cd "$input_dir" || exit

for bam in *_dedup.sorted.bam; do
  [ -e "$bam" ] || continue
  sample="${bam%%_dedup.sorted.bam}"

  echo "QC for $bam"
  $SAMTOOLS flagstat "$bam" > "${sample}.flagstat.txt"
  $SAMTOOLS idxstats "$bam" > "${sample}.idxstats.txt"
done

## MACS2 Peak Calling (NARROW)

echo "Running MACS2, narrow peaks"

export PYTHONPATH=/programs/macs2-2.2.9.1/lib64/python3.9/site-packages
export PATH=/programs/macs2-2.2.9.1/bin:$PATH

mkdir -p macs2_narrow
# If you have input files, add:
# input_bam="input_dedup.sorted.bam"
# and add: -c "$input_bam" \ below

for bam in *_dedup.sorted.bam; do
  [ -e "$bam" ] || continue
  base=$(basename "$bam" _dedup.sorted.bam)
  echo "Calling narrow peaks for $bam"
  macs2 callpeak \
    -t "$bam" \
    -f BAMPE \
    -g hs \
    -n "$base" \
    -q 0.01 \
    --nomodel \
    --keep-dup all \
    -B \
    --outdir macs2_narrow
done

## MACS2 Peak Calling (BROAD)

echo "Running MACS2, broad peaks"

mkdir -p macs2_broad
# If you have input files, add:
# input_bam="input_dedup.sorted.bam"
# and add: -c "$input_bam" \ below

for bam in *_dedup.sorted.bam; do
  [ -e "$bam" ] || continue
  base=$(basename "$bam" _dedup.sorted.bam)
  echo "Calling broad peaks for $bam"
  macs2 callpeak \
    -t "$bam" \
    -f BAMPE \
    -g hs \
    -n "$base" \
    -q 0.05 \
    --nomodel \
    --keep-dup all \
    --broad \
    -B \
    --outdir macs2_broad
done

## bedGraph -> BigWig (raw & FE)

echo "Converting bedGraph to BigWig and FE tracks"

for bdg in *_treat_pileup.bdg; do
  [ -e "$bdg" ] || continue
  base="${bdg%%_treat_pileup.bdg}"
  echo "Converting $bdg to BigWig"

  grep -v '^chrEBV' "$bdg" | sort -k1,1 -k2,2n > "${base}_treat_pileup_noEBVsorted.bdg"

  /workdir/sml327/HOPE/ChIP/2024_TFs/BACH2/figure/bedGraphToBigWig \
    "${base}_treat_pileup_noEBVsorted.bdg" \
    /workdir/sml327/HOPE/ChIP/hg38.chrom.sizes \
    "${base}_treat_pileup.bw"
done

for bdg in *_treat_pileup_noEBVsorted.bdg; do
  [ -e "$bdg" ] || continue
  base="${bdg%%_treat_pileup_noEBVsorted.bdg}"
  ctrl="${base}_control_lambda.bdg"

  echo "Generating FE BigWig for $base"

  grep -v '^chrEBV' "$ctrl" | sort -k1,1 -k2,2n > "${base}_control_lambda_noEBVsorted.bdg"

  macs2 bdgcmp -t "$bdg" -c "${base}_control_lambda_noEBVsorted.bdg" -m FE -o "${base}_FE.bdg"

  sort -k1,1 -k2,2n "${base}_FE.bdg" > "${base}_FE_sorted.bdg"
  /workdir/sml327/HOPE/ChIP/2024_TFs/BACH2/figure/bedGraphToBigWig \
    "${base}_FE_sorted.bdg" \
    /workdir/sml327/HOPE/ChIP/hg38.chrom.sizes \
    "${base}_FE.bw"
done

echo "ChIP-seq pipeline completed."