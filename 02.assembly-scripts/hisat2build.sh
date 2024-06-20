#!/bin/bash
#SBATCH --job-name=hisat2build
#SBATCH -N 1
#SBATCH --mem=80GB
#SBATCH --ntasks-per-node=20

EXPORT=/LUSTRE/apps/bioinformatica/hisat2/

export PATH=$PATH:$EXPORT

CPU=$SLURM_NPROCS

genome=$1 #

prefix=`basename ${genome%.f*}`

hisat2-build -p $CPU $genome $prefix

exit