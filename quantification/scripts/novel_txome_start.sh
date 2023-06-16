#!/bin/bash 
#SBATCH --job-name novel_tx_quant

wd=`pwd`
scriptdir="${wd}/scripts"
samples=(${wd}/data/raw/*fastq.gz)
simul_array_runs=50
run=$(echo 1-`date +%F`)

echo "`date` running salmon quant run $run"

mkdir -p log/${run}/${index_name}/

salmon_quant_jobid=()

salmon_quant_jobid+=($(sbatch --parsable \
  --mem=48G \
  --cpus-per-task=6 \
  --time=36:00:00 \
  --gres=tmpspace:50G \
  --array 1-${#samples[@]}%${simul_array_runs} \
  --job-name=${run}.quant_novel_txome \
  --output=log/${run}/%A_%a.out \
  ${scriptdir}/quant_novel_txome.sh \
))

echo "Novel salmon jobid: ${salmon_quant_jobid[@]}"
