#!/bin/bash

# Load modules
module load miniconda

mkdir -p ${outdir}/arcas-hla/json_files

 # Find all 'genotype.json' files and create symlinks in the destination directory
find "${outdir}/arcas-hla" -name '*genotype.json' -exec bash -c '
  outdir=$1
  shift
  for file; do
    # Generate symlink path
    symlink="${outdir}/arcas-hla/json_files/$(basename "$file")"

    # Check if a file with the same name already exists in the destination directory
    if [[ -e "$symlink" ]]; then
      echo "Warning: A file with the same name already exists: $symlink. Skipping..."
      continue
    fi

    # Create the symlink
    ln -s "$file" "$symlink"
    echo "Symlink created for $file -> $symlink"
  done
' bash "$outdir" {} +

arcasHLA merge --run "all" --i ${outdir}/arcas-hla/json_files --o ${outdir}/arcas-hla --verbose