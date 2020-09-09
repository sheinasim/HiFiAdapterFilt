#!/bin/bash

threads=1
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
    printf "No options specified. \nFiltering .bam files in working directoryi. \n"
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

echo "Converting .bam to .fastq on $(date)"

for x in `echo ${reads_pref}`
do
if [ ! -s ${outdir}/${x}.fastq ]
then
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
blastn -db $DBpath/pacbio_vectors_db -query ${outdir}/${x}.fasta -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue .01 -searchsp 1750000000000 -outfmt 6 > ${outdir}/${x}.contaminant.blastout &
wait
echo "Creating blocklist of reads to filter on $(date)."
cat ${outdir}/${x}.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' '{if ($3 >= 99) print $1}' | sort -u > ${outdir}/${x}.blocklist &
wait
echo "Removing adapter contaminated reads from .fastq on $(date)."
cat ${outdir}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" > ${outdir}/${x}.filt.fastq 
wait
echo "Finished on $(date)"
fi
done

