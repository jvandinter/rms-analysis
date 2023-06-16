#!/bin/bash

function check_annotation() {

    # Check whether to use a custom GTF / Rannot or use the prebuild
    # reference GTF / Rannot

    reference_annotation=$1
    reference_gtf=$2
    custom_annotation=$3
    custom_gtf=$4
    custom_annotation_package=$5
    reference_annotation_package=$6

    # Check to see whether custom GTF / Rannot was provided
    if [[ -f ${custom_gtf} ]]; then
        gtf=${custom_gtf}
    else
        gtf=${reference_gtf}
    fi

    if [[ -f ${custom_annotation} ]]; then
        rannot=${custom_annotation}
        annotation_package=${custom_annotation_package}
    else
        rannot=${reference_annotation}
        annotation_package=${reference_annotation_package}
    fi
        annot_name=$(basename ${rannot%.gtf_Rannot})
        annot_name=${annot_name#Homo_sapiens*.*}
}

function get_samples() {

    # Creates an array that contains all samples that will be analysed, including
    # additional arrays that contain the barcode and sample_ID

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
