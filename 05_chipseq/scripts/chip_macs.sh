#!/bin/bash

######################################################################
#
# Authors: J.T.vandinter-3@prinsesmaximacentrum.nl
# Authors: D.A.Hofman-3@prinsesmaximacentrum.nl
# Date: 03-06-2021
#
######################################################################

set -euo pipefail

# Load parameters from main script
threads=$((SLURM_CPUS_PER_TASK * 2))

module load macs/3.0.0

mkdir -p "${outdir}/macs3/"

cd ${TMPDIR}

macs3 --version

# # CTR-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/CTR-ChIP-H3K27ac" \
#         --name "CTR-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/755-CTR-ChIP-H3K27ac/755-CTR-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/756-CTR-ChIP-control/756-CTR-ChIP-control_dedup.bam"

# # NCI0075-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/NCI0075-ChIP-H3K27ac" \
#         --name "NCI0075-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/757-NCI0075-ChIP-H3K27ac/757-NCI0075-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/758-NCI0075-ChIP-control/758-NCI0075-ChIP-control_dedup.bam"

# # NCI0082-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/NCI0082-ChIP-H3K27ac" \
#         --name "NCI0082-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/759-NCI0082-ChIP-H3K27ac/759-NCI0082-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/760-NCI0082-ChIP-control/760-NCI0082-ChIP-control_dedup.bam"

# # NS129-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/NS129-ChIP-H3K27ac" \
#         --name "NS129-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/761-NS129-ChIP-H3K27ac/761-NS129-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/762-NS129-ChIP-control/762-NS129-ChIP-control_dedup.bam"

# # NS134-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/NS134-ChIP-H3K27ac" \
#         --name "NS134-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/763-NS134-ChIP-H3K27ac/763-NS134-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/764-NS134-ChIP-control/764-NS134-ChIP-control_dedup.bam"

# # RD-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RD-ChIP-H3K27ac" \
#         --name "RD-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/765-RD-ChIP-H3K27ac/765-RD-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/766-RD-ChIP-control/766-RD-ChIP-control_dedup.bam"

# # RH4-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RH4-ChIP-H3K27ac" \
#         --name "RH4-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/774-RH4-ChIP-H3K27ac/774-RH4-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# # RH4-ChIP-H3K27me3
# # broad
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RH4-ChIP-H3K27me3" \
#         --name "RH4-ChIP-H3K27me3" \
#         --qvalue 0.01 \
#         --treatment "${outdir}/picard/775-RH4-ChIP-H3K27me3/775-RH4-ChIP-H3K27me3_dedup.bam" \
#         --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# # RH4-ChIP-H3K36me3
# # broad
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RH4-ChIP-H3K36me3" \
#         --name "RH4-ChIP-H3K36me3" \
#         --qvalue 0.01 \
#         --treatment "${outdir}/picard/776-RH4-ChIP-H3K36me3/776-RH4-ChIP-H3K36me3_dedup.bam" \
#         --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# # RH4-ChIP-H3K4me1
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RH4-ChIP-H3K4me1" \
#         --name "RH4-ChIP-H3K4me1" \
#         --treatment "${outdir}/picard/777-RH4-ChIP-H3K4me1/777-RH4-ChIP-H3K4me1_dedup.bam" \
#         --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# # RH4-ChIP-H3K4me2
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RH4-ChIP-H3K4me2" \
#         --name "RH4-ChIP-H3K4me2" \
#         --treatment "${outdir}/picard/778-RH4-ChIP-H3K4me2/778-RH4-ChIP-H3K4me2_dedup.bam" \
#         --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# # RH4-ChIP-H3K4me3
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RH4-ChIP-H3K4me3" \
#         --name "RH4-ChIP-H3K4me3" \
#         --treatment "${outdir}/picard/779-RH4-ChIP-H3K4me3/779-RH4-ChIP-H3K4me3_dedup.bam" \
#         --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# # RH5-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RH5-ChIP-H3K27ac" \
#         --name "RH5-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/791-RH5-ChIP-H3K27ac/791-RH5-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/792-RH5-ChIP-control/792-RH5-ChIP-control_dedup.bam"

# # RMS008-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RMS008-ChIP-H3K27ac" \
#         --name "RMS008-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/793-RMS008-ChIP-H3K27ac/793-RMS008-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/794-RMS008-ChIP-control/794-RMS008-ChIP-control_dedup.bam"

# # RMS206-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RMS206-ChIP-H3K27ac" \
#         --name "RMS206-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/795-RMS206-ChIP-H3K27ac/795-RMS206-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/796-RMS206-ChIP-control/796-RMS206-ChIP-control_dedup.bam"

# # RMS209-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RMS209-ChIP-H3K27ac" \
#         --name "RMS209-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/797-RMS209-ChIP-H3K27ac/797-RMS209-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/798-RMS209-ChIP-control/798-RMS209-ChIP-control_dedup.bam"


# # RMS216-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RMS216-ChIP-H3K27ac" \
#         --name "RMS216-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/799-RMS216-ChIP-H3K27ac/799-RMS216-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/800-RMS216-ChIP-control/800-RMS216-ChIP-control_dedup.bam"

# # RMS237-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RMS237-ChIP-H3K27ac" \
#         --name "RMS237-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/801-RMS237-ChIP-H3K27ac/801-RMS237-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/802-RMS237-ChIP-control/802-RMS237-ChIP-control_dedup.bam"

# # RMS238-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RMS238-ChIP-H3K27ac" \
#         --name "RMS238-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/803-RMS238-ChIP-H3K27ac/803-RMS238-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/804-RMS238-ChIP-control/804-RMS238-ChIP-control_dedup.bam"

# # SCMC-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/SCMC-ChIP-H3K27ac" \
#         --name "SCMC-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/805-SCMC-ChIP-H3K27ac/805-SCMC-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/806-SCMC-ChIP-control/806-SCMC-ChIP-control_dedup.bam"

# # 7250-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/7250-ChIP-H3K27ac" \
#         --name "7250-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/827-7250-ChIP-H3K27ac/827-7250-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/825-7250-ChIP-control/825-7250-ChIP-control_dedup.bam"

# # 7250FP-ChIP-H3K27ac
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/7250FP-ChIP-H3K27ac" \
#         --name "7250FP-ChIP-H3K27ac" \
#         --treatment "${outdir}/picard/828-7250FP-ChIP-H3K27ac/828-7250FP-ChIP-H3K27ac_dedup.bam" \
#         --control "${outdir}/picard/826-7250FP-ChIP-control/826-7250FP-ChIP-control_dedup.bam"

# ## Narrow peaks

# # RH4-ChIP-BRD4
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RH4-ChIP-BRD4" \
#         --name "RH4-ChIP-BRD4" \
#         --treatment "${outdir}/picard/767-RH4-ChIP-BRD4/767-RH4-ChIP-BRD4_dedup.bam" "${outdir}/picard/768-RH4-ChIP-BRD4/768-RH4-ChIP-BRD4_dedup.bam" \
#         --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# # RH4-ChIP-CTCF
# macs3 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs3/RH4-ChIP-CTCF" \
#         --name "RH4-ChIP-CTCF" \
#         --treatment "${outdir}/picard/769-RH4-ChIP-CTCF/769-RH4-ChIP-CTCF_dedup.bam" "${outdir}/picard/770-RH4-ChIP-CTCF/770-RH4-ChIP-CTCF_dedup.bam" \
#         --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# RH4-ChIP-MED1
macs3 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs3/RH4-ChIP-MED1" \
        --name "RH4-ChIP-MED1" \
        --treatment "${outdir}/picard/782-RH4-ChIP-MED1/782-RH4-ChIP-MED1_dedup.bam" \
        --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# RH4-ChIP-MYCN
macs3 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs3/RH4-ChIP-MYCN" \
        --name "RH4-ChIP-MYCN" \
        --treatment "${outdir}/picard/783-RH4-ChIP-MYCN/783-RH4-ChIP-MYCN_dedup.bam" \
        --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# RH4-ChIP-MYOD1
macs3 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs3/RH4-ChIP-MYOD1" \
        --name "RH4-ChIP-MYOD1" \
        --treatment "${outdir}/picard/784-RH4-ChIP-MYOD1/784-RH4-ChIP-MYOD1_dedup.bam" \
        --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# RH4-ChIP-MYOG
macs3 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs3/RH4-ChIP-MYOG" \
        --name "RH4-ChIP-MYOG" \
        --treatment "${outdir}/picard/785-RH4-ChIP-MYOG/785-RH4-ChIP-MYOG_dedup.bam" \
        --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# RH4-ChIP-p300
macs3 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs3/RH4-ChIP-p300" \
        --name "RH4-ChIP-p300" \
        --treatment "${outdir}/picard/786-RH4-ChIP-p300/786-RH4-ChIP-p300_dedup.bam" \
        --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# RH4-ChIP-RAD21
macs3 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs3/RH4-ChIP-RAD21" \
        --name "RH4-ChIP-RAD21" \
        --treatment "${outdir}/picard/787-RH4-ChIP-RAD21/787-RH4-ChIP-RAD21_dedup.bam" \
        --control "${outdir}/picard/780-RH4-ChIP-control/780-RH4-ChIP-control_dedup.bam" "${outdir}/picard/781-RH4-ChIP-control/781-RH4-ChIP-control_dedup.bam"

# 7250-ChIP-BRD4
macs3 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs3/7250-ChIP-BRD4" \
        --name "7250-ChIP-BRD4" \
        --treatment "${outdir}/picard/823-7250-ChIP-BRD4/823-7250-ChIP-BRD4_dedup.bam" \
        --control "${outdir}/picard/825-7250-ChIP-control/825-7250-ChIP-control_dedup.bam"

# 7250FP-ChIP-BRD4
macs3 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs3/7250FP-ChIP-BRD4" \
        --name "7250FP-ChIP-BRD4" \
        --treatment "${outdir}/picard/824-7250FP-ChIP-BRD4/824-7250FP-ChIP-BRD4_dedup.bam" \
        --control "${outdir}/picard/826-7250FP-ChIP-control/826-7250FP-ChIP-control_dedup.bam"

# 7250-ChIP-PAX3-FOXO
macs3 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs3/7250-ChIP-PAX3-FOXO" \
        --name "7250-ChIP-PAX3-FOXO" \
        --treatment "${outdir}/picard/829-7250-ChIP-PAX3-FOXO1/829-7250-ChIP-PAX3-FOXO1_dedup.bam" \
        --control "${outdir}/picard/825-7250-ChIP-control/825-7250-ChIP-control_dedup.bam"

# 7250FP-ChIP-PAX3-FOXO1
macs3 \
  callpeak \
        --gsize "hs" \
        --format "BAM" \
        --outdir "${outdir}/macs3/7250FP-ChIP-PAX3-FOXO1" \
        --name "7250FP-ChIP-PAX3-FOXO1" \
        --treatment "${outdir}/picard/830-7250FP-ChIP-PAX3-FOXO1/830-7250FP-ChIP-PAX3-FOXO1_dedup.bam" \
        --control "${outdir}/picard/826-7250FP-ChIP-control/826-7250FP-ChIP-control_dedup.bam"





# # RD cell line
# apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/macs2-2.2.7.1.sif macs2 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs2" \
#         --name "RD_PAX3-FOXO1" \
#         --treatment "${outdir}/picard/SRR039129/SRR039129_dedup.bam" "${outdir}/picard/SRR039133/SRR039133_dedup.bam" \
#         --control "${outdir}/picard/SRR039130/SRR039130_dedup.bam"

# # RH cell line
# apptainer exec -B "/hpc:/hpc" --env LC_ALL=C.UTF-8 ${container_dir}/macs2-2.2.7.1.sif macs2 \
#   callpeak \
#         --gsize "hs" \
#         --format "BAM" \
#         --outdir "${outdir}/macs2" \
#         --name "RH4_PAX3-FOXO1" \
#         --keep-dup all \
#         --treatment "${outdir}/picard/SRR039132/SRR039132_dedup.bam" "${outdir}/picard/SRR039135/SRR039135_dedup.bam" \
#         --control "${outdir}/picard/SRR039131/SRR039131_dedup.bam" "${outdir}/picard/SRR039134/SRR039134_dedup.bam"
