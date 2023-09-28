#!/bin/bash 
#SBATCH --cpus-per-task 1
#SBATCH --mem=24G
#SBATCH --time=24:00:00

plus="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/02_riboseq/analysis/merged_psites/tumoroid/plus"
minus="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/02_riboseq/analysis/merged_psites/tumoroid/minus"
bedtools="/hpc/local/Rocky8/pmc_vanheesch/software/bedtools2-2.31.0/bin/bedtools"

## PLUS ##

$bedtools unionbedg -i ${plus}/*.sorted.bedgraph > ${plus}/RMS_merged.sorted.bedgraph
$bedtools unionbedg -i ${plus}/*FP*.sorted.bedgraph > ${plus}/RMS_FP_merged.sorted.bedgraph
$bedtools unionbedg -i ${plus}/*FN*.sorted.bedgraph > ${plus}/RMS_FN_merged.sorted.bedgraph

# RMS
awk '{ sum = 0; for (i = 4; i <= NF; i++) sum += $i; print $1"\t"$2"\t"$3"\t"sum/(NF-3)}' ${plus}/RMS_merged.sorted.bedgraph > ${plus}/RMS_merged_avg.bedgraph
awk '{sum=0; for (i=4; i<=NF; i++) sum+=$i; print $1"\t"$2"\t"$3"\t"sum}' ${plus}/RMS_merged.sorted.bedgraph > ${plus}/RMS_merged_summed.bedgraph

# FP ORG RMS
awk '{ sum = 0; for (i = 4; i <= NF; i++) sum += $i; print $1"\t"$2"\t"$3"\t"sum/(NF-3)}' ${plus}/RMS_FP_merged.sorted.bedgraph > ${plus}/RMS_FP_merged_avg.bedgraph
awk '{sum=0; for (i=4; i<=NF; i++) sum+=$i; print $1"\t"$2"\t"$3"\t"sum}' ${plus}/RMS_FP_merged.sorted.bedgraph > ${plus}/RMS_FP_merged_summed.bedgraph

# FN ORG RMS
awk '{ sum = 0; for (i = 4; i <= NF; i++) sum += $i; print $1"\t"$2"\t"$3"\t"sum/(NF-3)}' ${plus}/RMS_FN_merged.sorted.bedgraph > ${plus}/RMS_FN_merged_avg.bedgraph
awk '{sum=0; for (i=4; i<=NF; i++) sum+=$i; print $1"\t"$2"\t"$3"\t"sum}' ${plus}/RMS_FN_merged.sorted.bedgraph > ${plus}/RMS_FN_merged_summed.bedgraph

## MINUS ##

$bedtools unionbedg -i ${minus}/*.sorted.bedgraph > ${minus}/RMS_merged.sorted.bedgraph
$bedtools unionbedg -i ${minus}/*FP*.sorted.bedgraph > ${minus}/RMS_FP_merged.sorted.bedgraph
$bedtools unionbedg -i ${minus}/*FN*.sorted.bedgraph > ${minus}/RMS_FN_merged.sorted.bedgraph

# RMS
awk '{ sum = 0; for (i = 4; i <= NF; i++) sum += $i; print $1"\t"$2"\t"$3"\t"sum/(NF-3)}' ${minus}/RMS_merged.sorted.bedgraph > ${minus}/RMS_merged_avg.bedgraph
awk '{sum=0; for (i=4; i<=NF; i++) sum+=$i; print $1"\t"$2"\t"$3"\t"sum}' ${minus}/RMS_merged.sorted.bedgraph > ${minus}/RMS_merged_summed.bedgraph

# FP ORG RMS
awk '{ sum = 0; for (i = 4; i <= NF; i++) sum += $i; print $1"\t"$2"\t"$3"\t"sum/(NF-3)}' ${minus}/RMS_FP_merged.sorted.bedgraph > ${minus}/RMS_FP_merged_avg.bedgraph
awk '{sum=0; for (i=4; i<=NF; i++) sum+=$i; print $1"\t"$2"\t"$3"\t"sum}' ${minus}/RMS_FP_merged.sorted.bedgraph > ${minus}/RMS_FP_merged_summed.bedgraph

# FN ORG RMS
awk '{ sum = 0; for (i = 4; i <= NF; i++) sum += $i; print $1"\t"$2"\t"$3"\t"sum/(NF-3)}' ${minus}/RMS_FN_merged.sorted.bedgraph > ${minus}/RMS_FN_merged_avg.bedgraph
awk '{sum=0; for (i=4; i<=NF; i++) sum+=$i; print $1"\t"$2"\t"$3"\t"sum}' ${minus}/RMS_FN_merged.sorted.bedgraph > ${minus}/RMS_FN_merged_summed.bedgraph
