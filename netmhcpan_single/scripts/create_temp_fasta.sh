#!/bin/bash

outdir="/hpc/pmc_vanheesch/projects/Jip/rms_analysis/netmhcpan_single/analysis"
processed_fasta="${outdir}/netMHCpan/RMS_neoantigen_strict/peptides.fasta"
temp_fasta="${outdir}/netMHCpan/RMS_neoantigen_strict/temp"
name_table="${outdir}/netMHCpan/RMS_neoantigen_strict/temp_names.tsv"

# Create a temp name file and fasta file because NetMHCPan cant handle long input names in the fasta
counter=1

while IFS= read -r line; do
    if [[ $line =~ ^\>(.*)$ ]]; then
        echo $counter
        name="${BASH_REMATCH[1]}"
        new_name=$(printf "%015d" $counter)
        echo -e "$name\t$new_name" >> "$name_table"
        if [[ $name =~ (.*)__8__(.*) ]]; then
            echo ">$new_name" >> "${temp_fasta}_8.fasta"
            ((counter++))
        elif [[ $name =~ (.*)__9__(.*) ]]; then
            echo ">$new_name" >> "${temp_fasta}_9.fasta"
            ((counter++))
        elif [[ $name =~ (.*)__10__(.*) ]]; then
            echo ">$new_name" >> "${temp_fasta}_10.fasta"
            ((counter++))
        elif [[ $name =~ (.*)__11__(.*) ]]; then
            echo ">$new_name" >> "${temp_fasta}_11.fasta"
            ((counter++))
        elif [[ $name =~ (.*)__12__(.*) ]]; then
            echo ">$new_name" >> "${temp_fasta}_12.fasta"
            ((counter++))
        fi
    else 
        if [[ $name =~ (.*)__8__(.*) ]]; then
            echo "$line" >> "${temp_fasta}_8.fasta"
        elif [[ $name =~ (.*)__9__(.*) ]]; then
            echo "$line" >> "${temp_fasta}_9.fasta"
        elif [[ $name =~ (.*)__10__(.*) ]]; then
            echo "$line" >> "${temp_fasta}_10.fasta"
        elif [[ $name =~ (.*)__11__(.*) ]]; then
            echo "$line" >> "${temp_fasta}_11.fasta"
        elif [[ $name =~ (.*)__12__(.*) ]]; then
            echo "$line" >> "${temp_fasta}_12.fasta"
        fi
    fi
done < "${processed_fasta}"

tail "${name_table}"
