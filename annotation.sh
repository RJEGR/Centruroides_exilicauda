#!/bin/sh
## Directivas
#SBATCH --job-name=Trinotate
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

FILE=$1 #rnaspades.fasta
# prefix=`basename ${FILE%.f*}`

EXPORT=/LUSTRE/apps/bioinformatica/TransDecoder-v5.7.0/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/eggnog-mapper-master
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/tmhmm-2.0c/bin
export PATH=$PATH:$EXPORT

# Load depedencies for trinotate (such as blast, hmmpress, etc)
module load trinotate

TransDecoder.LongOrfs -t $FILE
TransDecoder.Predict -t $FILE --cpu $SLURM_NPROCS


#* to run ALL supported computes
TRINOTATE_DATA_DIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Trinotate/TRINOTATE_DB/

Trinotate --db sqlite.db --run ALL --CPU $SLURM_NPROCS --transcript_fasta $FILE --transdecoder_pep ${FILE}.transdecoder.pep --trinotate_data_dir  $TRINOTATE_DATA_DIR --use_diamond

exit