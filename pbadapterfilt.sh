#!/bin/bash

threads=1

unset name

while getopts 'b:t:' option
do
case "${option}"
in
b) bamprefix=${OPTARG};;
t) threads=${OPTARG} ;;
?) echo "Usage: $0 [ -b .bam prefix ] [ -t Number of threads for blastn. Default=1 ]"
esac
done

if ((OPTIND == 1))
then
    printf "No options specified. \nFiltering .bam files in working directory."
fi

shift $((OPTIND - 1))

reads=$(ls ${bamprefix}*bam | sed 's/\.bam//')

for x in `echo $reads`
do
if [ ! -s ${x}.fasta ]
then    
bamtools convert -format fasta -in ${x}.bam -out ${x}.fasta
fi
done

for x in `echo $reads`
do
if [ ! -s ${x}.fastq ]
then
bamtools convert -format fastq -in ${x}.bam -out ${x}.fastq
fi
done

for x in `echo $reads`
do
if [ ! -s ${x}.contaminant.blastout ]
then
blastn -db DB/pacbio_vectors_db -query ${x}.fasta -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue .01 -searchsp 1750000000000 -outfmt 6 >${x}.contaminant.blastout
fi
done

for x in `echo $reads`
do
if [ ! -s ${x}.blocklist ]
then
cat ${x}.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' '{if ($3 > 95) print $1}' | sort -u >${x}.blocklist
fi
done

for x in `echo $reads`
do
if [ ! -s ${x}.filt.fastq ]
then
cat ${x}.fastq | paste - - - - | grep -v -f ${x}.blocklist -F | tr "\t" "\n" > ${x}.filt.fastq
fi
done


