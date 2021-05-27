#!/bin/bash

threads=8
outdir=$(pwd)
DBpath=$(echo $PATH | sed 's/:/\n/g' | grep "HiFiAdapterFilt/DB" | head -n 1)

unset name

while getopts 'b:t:o:h' option
do
case "${option}"
in
b) bamprefix=${OPTARG};;
t) threads=${OPTARG} ;;
o) outdir=${OPTARG} ;;
h) echo "Usage: $0 [ -b .bam prefix ] [ -t Number of threads for blastn. Default=1 ] [ -o outdirectory prefix Default=. ]" ;;
?) echo "Usage: $0 [ -b .bam prefix ] [ -t Number of threads for blastn. Default=1 ] [ -o outdirectory prefix Default=. ]"
esac
done

if ((OPTIND == 1))
then
    printf "No options specified. \nFiltering files in working directory. \n"
fi

shift $((OPTIND - 1))

## Set variables based on options

reads_rp=$(ls ${bamprefix}*bam | sed 's/\.bam//')
reads_pref=$(ls ${bamprefix}*bam | sed 's/\.bam//' | rev | cut -d'/' -f 1 | rev)
read_path=$(dirname ${bamprefix}*.bam)
read_path_str=$(echo ${read_path} | cut -d" " -f 1)

## Create out directory if necessary

if [ ! -d ${outdir} ]
then 
	mkdir ${outdir}
fi

## Convert .bam to .fastq

for x in `echo ${reads_pref}`
do
if [ ! -s ${outdir}/${x}.fastq ]
then
    echo "Converting .bam to .fastq on $(date)"
    bamtools convert -format fastq -in ${read_path_str}/${x}.bam -out ${outdir}/${x}.fastq &
fi
done

## Convert .bam to .fasta

for x in `echo ${reads_pref}`
do
if [ ! -s ${outdir}/${x}.fasta ]
then
	echo "Converting .bam to .fasta on $(date)."    
	bamtools convert -format fasta -in ${read_path_str}/${x}.bam -out ${outdir}/${x}.fasta &
	wait
	echo "Identifying reads with adapter contamination on $(date)."
	blastn -db $DBpath/pacbio_vectors_db -query ${outdir}/${x}.fasta -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue .01 -searchsp 1750000000000 -outfmt 6 > ${outdir}/${x}.contaminant.blastout &
	wait
	echo "Creating blocklist of reads to filter on $(date)."
	cat ${outdir}/${x}.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' '{if (($2 ~ /NGB00972/ && $3 >= 97 && $4 >= 44) || ($2 ~ /NGB00973/ && $3 >= 97 && $4 >= 34)) print $1}' | sort -u > ${outdir}/${x}.blocklist &  
	wait
	echo "Removing adapter contaminated reads from .fastq on $(date)."
	cat ${outdir}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" > ${outdir}/${x}.filt.fastq 
else
	echo "Identifying reads with adapter contamination on $(date)."
	blastn -db $DBpath/pacbio_vectors_db -query ${outdir}/${x}.fasta -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue .01 -searchsp 1750000000000 -outfmt 6 > ${outdir}/${x}.contaminant.blastout &
	wait
	echo "Creating blocklist of reads to filter on $(date)."
	cat ${outdir}/${x}.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' '{if (($2 ~ /NGB00972/ && $3 >= 97 && $4 >= 44) || ($2 ~ /NGB00973/ && $3 >= 97 && $4 >= 34)) print $1}' | sort -u > ${outdir}/${x}.blocklist &
	wait
	echo "Removing adapter contaminated reads from .fastq on $(date)."
	cat ${outdir}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" > ${outdir}/${x}.filt.fastq
	wait
fi
done

f=`cat ${outdir}/${x}.blocklist | wc -l` #number of adapter contaminated reads
r1=`cat ${outdir}/${x}.fastq | wc -l` 
r2=`awk -v r1=$r1 'BEGIN{ans=r1/4; print ans}'` #number of ccs reads
p1=`awk -v n1=$r2 -v n2=$f 'BEGIN{ans=n2/n1*100; print ans}'` #proportion of adapter contaminated reads
r3=`awk -v r2=$r2 -v f=$f 'BEGIN{ans=r2-f; print ans}'` #number of reads retained
p2=`awk -v p1=$p1 'BEGIN{ans=1-p1; print ans}'` #proportion of reads retained

touch ${outdir}/${x}.stats

echo "For the" ${x} "dataset:" >>${outdir}/${x}.stats 
echo "" >>${outdir}/${x}.stats
echo "Number of ccs reads:" $r2 >>${outdir}/${x}.stats
echo "Number of adapter contaminated ccs reads:" $f "("$p1"% of total)" >>${outdir}/${x}.stats
echo "Number of ccs reads retained:" $r3 "("$p2"% of total)" >>${outdir}/${x}.stats
echo "" >>${outdir}/${x}.stats
echo "Finished on $(date)" >>${outdir}/${x}.stats
