#!/bin/sh
## Directivas
#SBATCH --job-name=rsem
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

thread_count=$SLURM_NPROCS

# 2. Prepare reference for rsem

ref=$1 # reference.fasta file

rsem_prefix=idx # conveniently as $ref_index_name

rsem-prepare-reference $ref $rsem_prefix

# 2.1) Convert bam to rsem 

for i in $(ls *.sorted.bam);
do
bam_for_rsem=${i%.bam}.rsem
convert-sam-for-rsem -p $thread_count $i $bam_for_rsem
done

# 3) Estimate abundances
# RSEM Parameters

fragment_length=200
fragment_std=80
fraglength_info_txt="--fragment-length-mean $fragment_length --fragment-length-sd $fragment_std"

paired_flag_text="--paired-end"  
no_qualities_string=""
keep_intermediate_files_opt="--keep-intermediate-files"

SS_opt="--forward-prob 1.0"
rsem_bam_flag="--no-bam-output"

# Estimate abundance:

for i in $(ls *.rsem.bam);
do
output_prefix=${i%.sorted.rsem.bam}
rsem-calculate-expression $no_qualities_string $paired_flag_text -p $thread_count $fraglength_info_txt $keep_intermediate_files_opt $SS_opt $rsem_bam_flag --bam $i $rsem_prefix $output_prefix

done

exit
