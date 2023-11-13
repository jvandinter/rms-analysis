#!/usr/bin/env python3

"""Prepare the input to run MHCFlurry and remove non-tumor-specific peptides.

Authors: l.w.kok-15@prinsesmaximacentrum.nl
Date: 05-11-2021

This script reads proteins from an input .xls or .fasta file and extracts all 
peptides with a length from 8-11 amino acids from the file. Then it removes 
peptides occuring in the human proteome. The remaining (tumor-specific)
peptides are written to a .csv output file formatted such that it can be used as 
input by MHCFlurry 

Input arguments:
    argv[1] - the input .xls or .fasta file. The .xls file should contain the 
              columns 'Protein' and 'ORF_id'.
    argv[2] - file to which the output, to be used by MHCFlurry, should be
              written.
    argv[3] - directory with either .faa and .fasta files containing the human
              proteome, or the proteome_kmers files created earlier by this
              script.
    argv[4] - the HLA allele for which binding should be predicted - 
              depracated
"""

from sys import argv
import pandas as pd
from pathlib import Path
import json


def get_proteome_kmers(prot_dir, len_lim):
    """Get kmers from protein sequences in the proteome.

    Keyword arguments:
        prot_dir -- path, directory with the proteome in fasta files, and/or
            proteome kmer files created earlier by this function.
        len_lim -- list of two ints, minimum and maximum length of the kmers.

    Returns:
        kmers_ref -- list of sets, each set contains kmers of a specific 
            length. The list is sorted by the kmer length in the sets.

    Info:
        All produced kmers are saved per length in prot_dir as txt files.
    """     
    prot_lim = []
    kmers_ref = []
    for l in range(len_lim[0], len_lim[1]+1):
        if not (prot_dir / "proteome_kmers_{}.txt".format(l)).is_file():
            prot_lim.append(l)
    if prot_lim:
        proteome = set()
        prot_files = [f for p in ["*.faa", "*.fa", "*.fasta"] for f in prot_dir.glob(p) 
            if not f.stem.startswith(".")]
        prot_kmers = list(zip(*[get_k_mers_per_len(s[0], prot_lim[0], 
            prot_lim[-1]) for f in prot_files for s in read_fasta(f)]))
    for i in range(len_lim[0], len_lim[1]+1):
        if not (prot_dir / "proteome_kmers_{}.txt".format(i)).is_file():
            k_mers = set().union(*prot_kmers[i-prot_lim[0]])
            kmers_ref.append(k_mers)
            file_name = "proteome_kmers_{}.txt".format(i)
            with open(prot_dir / file_name, "w+") as f_open:
                json.dump(list(k_mers), f_open)
        else:
            with open(prot_dir / "proteome_kmers_{}.txt".format(i)) as k_open:
                kmers_ref.append(set(json.load(k_open)))
    return kmers_ref    


def get_k_mers_per_len(seq, min_len, max_len):
    """Finds all k-mers in given length range in given string.

    Input parameters:
        seq -- string, sequence from which the k-mers have to be extracted
        min_len -- int, minimum length of the k-mers
        max_len -- int, maximum length of the k-mers

    Returns:
        k_mers -- list of sets, each set contains k_mers of one length
    """
    k_mers = [set() for x in range(min_len, max_len + 1)]
    for i in range(len(seq)):
        for j in range(i + min_len, min(i+max_len, len(seq))+1):
            k_mers[j-min_len-i].add(seq[i:j])
    return k_mers


def get_k_mers(seq, seq_id, min_len, max_len, ref_kmer):
    """Finds k-mers of given lengths that do not occur in the reference.

    Input parameters:
        seq -- string, sequence from which the k-mers have to be extracted. Only 
            kmers not occuring in ref_kmers will be kept.
        seq_id -- string, id of the sequence.
        min_len -- int, minimum length of the k-mers.
        max_len -- int, maximum length of the k-mers.
        ref_kmer -- list of sets, each set contains kmers of a specific 
            length. The list is sorted by the kmer length in the sets.

    Returns:
        k_mers -- list of tuples, tuples are the k_mer, n-flanking region, 
            c-flanking region and a seq_id derived unique id.
    """
    count = 0
    k_mers = []
    len_seq = len(seq)
    for i in range(len(seq)):
        for j in range(i + min_len, min(i+max_len, len(seq))+1):
            if seq[i:j] not in ref_kmer[j-min_len-i]:
                k_mer_id = "{}__{}__{}".format(seq_id, j-i, count)
                min_n_flank = max(0, i-15)
                max_c_flank = min(len_seq, j+15)
                k_mers.append((seq[i:j], seq[min_n_flank:i], 
                    seq[j:max_c_flank], k_mer_id))
        count += 1
    return k_mers


def read_fasta(fasta_file):
    """Yields sequences and sequence id from a fasta file

    Input parameters:
        fasta_file -- path, path to the input fasta file.

    Yields:
        sequence -- string
        sequence_id -- string, the id of the sequence
    """
    sequence = ""
    with open(fasta_file) as f_open:
        for line in f_open:
            if not line.strip():
                continue
            if line.startswith(">"):
                if sequence:
                    yield sequence, sequence_id
                    sequence = ""
                sequence_id = line[1:].strip()
            else:
                sequence += line.strip()
    yield sequence, sequence_id


def process_protein(seq, seq_id, len_lim, f_out, ref_kmers, allele):
    """Extract k_mers from a sequence and save them in csv format.

    Keyword arguments:
        seq -- string, the sequence from which k_mers will be extracted. Only
            k_mers not occuring in ref_kmers will be kept.
        seq_id -- string, id of the sequence.
        len_lim -- list of two ints, represent the min and max k_mer length.
        f_out -- open file, file to save the k_mers. 
        ref_kmers -- list of sets, each set contains kmers of a specific 
            length. The list is sorted by the kmer length in the sets.
        allele -- string, MHC-I allele for which predictions should be made.
    """
    k_mers = get_k_mers(seq, seq_id, *len_lim, ref_kmers)
    [f_out.write(">{}\n{}\n".format(k_mer[3], k_mer[0])) for k_mer in k_mers]


if __name__ == "__main__":
    pep_lim = [8, 12]
    in_file = Path(argv[1])
    prot_kmers = get_proteome_kmers(Path(argv[3]), pep_lim)
    hla_allele = ""
    with open(Path(argv[2]), "w+") as out_open:
        if in_file.suffix.startswith(".xls"):
            proteins = pd.read_excel(argv[1], engine='openpyxl')
            proteins.apply(lambda x: process_protein(x['Protein'], x['ORF_id'], 
            pep_lim, out_open, prot_kmers, hla_allele), axis=1)
        elif in_file.suffix == ".fasta":
            for pro, pro_id in read_fasta(in_file):
                process_protein(pro, pro_id, pep_lim, out_open, prot_kmers, 
                    hla_allele)
        else:
            raise IOError("{} is not an accepted file suffix".format(
                in_file.suffix))
