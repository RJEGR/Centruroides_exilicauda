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

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/tmhmm-2.0c/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/Trinotate/util/trinotateSeqLoader/
export PATH=$PATH:$EXPORT

# Load depedencies for trinotate (such as blast, hmmpress, etc)
module load trinotate

# Load this for run mapper.pl

module load conda-2024
source activate base
conda activate trinotate2

TransDecoder.LongOrfs -t $FILE
TransDecoder.Predict -t $FILE --cpu $SLURM_NPROCS


#* to run ALL supported computes
TRINOTATE_DATA_DIR=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Trinotate/TRINOTATE_DB/

rnaspades_trans_map

grep "^>" rnaspades.fasta | awk '{gsub(/_i[0-9]|>/,""); print}' > genes
grep "^>" rnaspades.fasta | awk '{gsub(/>/,""); print}' > trans

paste -d"\t" genes trans >  gen_tr_map

Trinotate --db sqlite.db --init --gene_trans_map gen_tr_map --transcript_fasta rnaspades.fasta --transdecoder_pep rnaspades.fasta.transdecoder.pep


Trinotate --db sqlite.db --run ALL --CPU $SLURM_NPROCS --transcript_fasta $FILE --transdecoder_pep ${FILE}.transdecoder.pep --trinotate_data_dir  $TRINOTATE_DATA_DIR --use_diamond

exit

module load conda-2024
source activate base
conda activate trinotate2
# conda deactivate
# conda activate trinotate_test

emapper.py -i rnaspades.fasta.transdecoder.pep --cpu 20 --data_dir /LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Trinotate/TRINOTATE_DB//EGGNOG_DATA_DIR -o eggnog_mapper --override



module load trinotate

Trinotate --db sqlite.db --LOAD_EggnogMapper eggnog_mapper.emapper.annotations

Trinotate --db sqlite.db --LOAD_custom_blast anchnoDB.diamond.blastp.outfmt6 --blast_type blastp --custom_db_name anchnoDB

Trinotate --db sqlite.db --report --incl_pep > Trinotate.xls


/LUSTRE/apps/bioinformatica/Trinotate/util/Trinotate_report_writer.pl --sqlite sqlite.db

# 
EXPORT=/LUSTRE/apps/bioinformatica/Trinotate/util/
export PATH=$PATH:$EXPORT

extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate.xls --gene > gene_ontology
extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate.xls --trans > trans_ontology

# if emtpy result, run init
# same issue https://www.biostars.org/p/9598219/

# check loading checkpoints in __trinotate_run_checkpts/, if neeeded run:
#Trinotate --db <sqlite.db> --LOAD_swissprot_blastp <file.outfmt6>
# Trinotate --db <sqlite.db> --LOAD_pfam <file>
#Trinotate --db <sqlite.db> --LOAD_signalp6 <file>
#Trinotate --db <sqlite.db> --LOAD_EggnogMapper <file>
# Trinotate --db <sqlite.db> --LOAD_tmhmmv2 <file> 
# Trinotate --db <sqlite.db> --LOAD_deeptmhmm <file.gff3>