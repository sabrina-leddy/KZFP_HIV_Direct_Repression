#!/bin/bash

# Exit on error
set -e

# Pipeline for strand-aware bigWigs

## Configuration
# deepTools
export PATH="/programs/deeptools-3.5.5/bin:$PATH"
export PYTHONPATH="/programs/deeptools-3.5.5/lib64/python3.9/site-packages:/programs/deeptools-3.5.5/lib/python3.9/site-packages"
# Samtools
SAMTOOLS="/programs/samtools-1.9-r9/bin/samtools"
# TMPDIR
mkdir -p "/workdir/sml327/tmp"
export TMPDIR="/workdir/sml327/tmp"

# Parameters
THREADS=4
BINSIZE=10
NORMALIZATION="CPM"

## bigWig prep

echo "Generating strand-aware BigWigs"

BAMS=(*.bam)

if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "ERROR: No .bam files found in this directory."
    exit 1
fi

# Forward strand
parallel -j 2 "
    base=\$(basename {} .bam);
    echo 'Forward: \$base';
    bamCoverage \
        -b {} \
        -o \${base}_forward_${NORMALIZATION}.bw \
        --filterRNAstrand forward \
        --normalizeUsing ${NORMALIZATION} \
        --binSize ${BINSIZE} \
        -p ${THREADS}
" ::: "${BAMS[@]}"

# Reverse strand
parallel -j 2 "
    base=\$(basename {} .bam);
    echo 'Reverse: \$base';
    bamCoverage \
        -b {} \
        -o \${base}_reverse_${NORMALIZATION}.bw \
        --filterRNAstrand reverse \
        --normalizeUsing ${NORMALIZATION} \
        --binSize ${BINSIZE} \
        -p ${THREADS}
" ::: "${BAMS[@]}"

echo "BigWig creation completed."