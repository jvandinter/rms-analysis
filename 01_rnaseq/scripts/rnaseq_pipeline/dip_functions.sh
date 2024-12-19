#!/bin/bash

function info() {
    echo "INFO: $@" >&2
}
function error() {
    echo "ERR:  $@" >&2
}
function fatal() {
    echo "ERR:  $@" >&2
    exit 1
}

function get_samples() {

    wd=$1
    data_folder=$2

    # Check whether the files are in correct format
    if [[ $(ls ${wd}/data/raw/*.fastq.gz | wc -l) -eq 0 ]]; then

        echo "ERROR: Pipeline only accepts .fastq.gz. Please change fq extension to fastq."
        exit

    else
        fastq_names=($(find "${wd}/data/raw" -maxdepth 1 -name "*_R1*" -exec basename {} \; 2> >(grep -v 'Permission denied' >&2) | sort -u ))

        for i in ${!fastq_names[@]}; do

            fastq_files[i]=$(find "${data_folder}" -name ${fastq_names[i]} 2> >(grep -v 'Permission denied' >&2))

        done

    fi

# Initiate arrays
    sample_ids=()
    barcode_ids=()
    samples=()

# Get sample IDs and barcodes from fastq files
    for i in ${!fastq_files[@]}; do

        sample_ids[i]=$(basename ${fastq_files[i]} | cut -f 1 -d "_")
        barcode_ids[i]=$(basename ${fastq_files[i]} | cut -f 2 -d "_")
        samples[i]=$(basename ${fastq_files[i]})

    done
}