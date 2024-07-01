#!/bin/bash
#SBATCH --job-name=Trinity
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

THREADS=$SLURM_NPROCS
MEM=${SLURM_MEM_PER_NODE}G
DIR=$SLURM_SUBMIT_DIR

module load trinityrnaseq-v2.15.1

# Ejemplo
# Trinity --seqType fq --max_memory 100G --left S1.F.fq.gz,S2.F.fq.gz,S3.F.fq.gz --right S1.R.fq.gz,S2.R.fq.gz,S3.R.fq.gz --CPU 20

,
,

,



,

--

rnaspades.py --pe1-1 CSI_M_PR_AD_24_1_F.P.qtrim.fq.gz --pe1-2 CSI_M_PR_AD_24_1_R.P.qtrim.fq.gz \
             --pe1-1 CSI_M_PR_AD_24_2_F.P.qtrim.fq.gz --pe1-2 CSI_M_PR_AD_24_2_R.P.qtrim.fq.gz \

             --pe1-1 GLO_M_PO_AD_24_P_F.fastq.gz --pe1-2 GLO_M_PO_AD_24_P_R.fastq.gz \
             --pe1-1 GLO_M_PR_AD_24_P_F.fastq.gz --pe1-2 GLO_M_PR_AD_24_P_R.fastq.gz \
             --pe1-1 LOP_M_PO_AD_24_P_F.fastq.gz --pe1-2 LOP_M_PO_AD_24_P_R.fastq.gz \
             --pe1-1 LOP_M_PR_AD_24_P_F.fastq.gz --pe1-2 LOP_M_PR_AD_24_P_R.fastq.gz \
             --pe1-fr -t $THREADS -m $MEM --cov-cutoff auto -o $DIR/rnaspades
