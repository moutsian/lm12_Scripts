#!/bin/bash
chrom=$1;
bsub -q normal -J qu${chrom} -R"select[mem>10000] rusage[mem=10000]" -M10000 -o quant${chrom}.out -e quant${chrom}.err /software/R-3.3.0/bin/Rscript ~lm12/Scripts/Rscript_for_aligned_read_pairs.R "${chrom}"

