#ssh rgomez@omica.cicese.mx 

#awk 'BEGIN{OFS=","} {$1=$1; print}' samples.txt | cut -f3 -d","

fasta=$1
samples=$2
#left=$2
#right=$3
out=${fasta%.fa}

 #for f in $(awk 'BEGIN{OFS=","} {$1=$1; print}' samples.txt | cut -f3 -d","); 
 #do ln -s $f .; done

 #for f in $(awk 'BEGIN{OFS=","} {$1=$1; print}' samples.txt | cut -f4 -d","); 
 #do ln -s $f .; done

 #cat *_1.P.qtrim.fq.gz > tmp.1.P.qtrim.fq.gz
 #cat *_2.P.qtrim.fq.gz > tmp.2.P.qtrim.fq.gz

for f in $(awk 'BEGIN{OFS=","} {$1=$1; print}' $samples | cut -f3 -d","); do cat $f >> tmp.1.P.qtrim.fq.gz .; done

for f in $(awk 'BEGIN{OFS=","} {$1=$1; print}' $samples | cut -f4 -d","); do cat $f >> tmp.2.P.qtrim.fq.gz .; done

left="ls tmp.1.*"
right="ls tmp.2.*

transrate --assembly $fasta --left $left --right $right --threads $SLURM_NPROCS --output ${out}_transrate

exit