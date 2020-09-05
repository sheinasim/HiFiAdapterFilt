#!/bin/bash

threads=1
outdir=$(pwd)

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

## Convert .bam to .fastq

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
bamtools convert -format fasta -in ${read_path_str}/${x}.bam -out ${outdir}/${x}.fasta
fi
done

## Blast raw reads to PB adapter database

for x in `echo ${reads_pref}`
do
if [ ! -s ${outdir}/${x}.contaminant.blastout ]
then
blastn -db $BLASTDB/pacbio_vectors_db -query ${outdir}/${x}.fasta -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue .01 -searchsp 1750000000000 -outfmt 6 > ${outdir}/${x}.contaminant.blastout
fi
done

## Create blocklist of reads to remove

for x in `echo ${reads_pref}`
do
if [ ! -s ${outdir}/${x}.blocklist ]
then
cat ${outdir}/${x}.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' '{if ($3 > 95) print $1}' | sort -u > ${outdir}/${x}.blocklist
fi
done

## Filter adapter contaminated reads from .fastq

for x in `echo ${reads_pref}`
do
if [ ! -s ${outdir}/${x}.filt.fastq ]
then
cat ${outdir}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" > ${outdir}/${x}.filt.fastq
fi
done


