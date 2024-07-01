#!/bin/sh
## Directivas
#SBATCH --job-name=fcount
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

export=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/subread-2.0.6-Linux-x86_64/bin
export PATH=$PATH:$EXPORT

featureCounts -T $SLURM_NPROCS -a rnaspades.fasta.transdecoder.gtf -o feature_counts.txt *bam

exit

ls -1 *_transcripts.gtf > stringtie_gtf_list.txt

stringtie --rf --merge -p 20 -o transcripts.gtf stringtie_gtf_list.txt



