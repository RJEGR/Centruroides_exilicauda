#!/bin/sh
## Directivas
#SBATCH --job-name=bbuild
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

EXPORT=/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:$EXPORT

# 1) INDEX THE REFERENCE

reference=$1 # reference.fasta
ref_index_name=idx # ${reference%.cdna}

# Build Bowtie2 index if not already present
if [ ! -f "$ref_index_name.1.bt2" ]; then
    bowtie2-build --threads $SLURM_NPROCS $reference $ref_index_name
else
    echo "index '$ref_index_name' already exists."
fi

exit