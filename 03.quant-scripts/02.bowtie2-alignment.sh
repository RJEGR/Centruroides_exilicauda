#!/bin/sh
## Directivas
#SBATCH --job-name=bowtie2
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

EXPORT=/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT


reference=$1 # reference.fasta
FirstPair=$2 # File_1.P.qtrim.fq.gz, the script will detect File_2.P.qtrim.fq.gz
ref_index_name=idx # conveniently as $ref_index_name

thread_count=$SLURM_NPROCS

# Parameters for alignment

aligner_params="--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200"
read_type="-q"
max_ins_size=800

base="${FirstPair%*_1.P.qtrim.fq.gz}"
left_file=${base}_1.P.qtrim.fq.gz
right_file=${base}_2.P.qtrim.fq.gz

bam_file=${base}.sorted.bam
met_file=${base}.met.txt
met_stderr=${base}.stderr # too memory  --met-stderr $met_stderr


# 1) Generate index from reference

# Build Bowtie2 index if not already present
if [ ! -f "$bam_file" ]; then
    
    echo "Aligning reads back to reference"

    bowtie2 --met-file $met_file $aligner_params $read_type -X $max_ins_size \
    -x $ref_index_name -1 $left_file -2 $right_file -p $thread_count 2> ${base}.stderr | \
    samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o $bam_file
else
    echo "File $bam_file already exists."

fi

exit