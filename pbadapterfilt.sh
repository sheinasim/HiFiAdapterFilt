#!/bin/bash

version=2.0.1
prefix=all
threads=8
outdir=$(pwd)
DBpath=$(echo $PATH | sed 's/:/\n/g' | grep "HiFiAdapterFilt/DB" | head -n 1)
adapterlength=44
pctmatch=97

USAGE="Usage: $0 [ -p sequence file Prefix ] [ -l minimum match Length to filter. Default=44 ] [ -m minimum Match percentage to filter. Default=97]  [ -t number of Threads     for blastn. Default=8 ] [ -o Outdirectory prefix Default=. ]"

unset name


while getopts ':p:l:m:t:o:cvh-' OPT
do
  if [ "$OPT" = "-" ]; then   # long option: reformulate OPT and OPTARG
    OPT="${OPTARG%%=*}"       # extract long option name
    OPTARG="${OPTARG#$OPT}"   # extract long option argument (may be empty)
    OPTARG="${OPTARG#=}"      # if long option argument, remove assigning `=`
  fi
case "${OPT}"
in
p | prefix ) prefix=${OPTARG};;
l | adapterlength ) adapterlength=${OPTARG};;
m | pctmatch ) pctmatch=${OPTARG};;
t | threads ) threads=${OPTARG} ;;
o | outdir ) outdir=${OPTARG} ;;
c | cite ) echo "Sim, S.B., Corpuz, R.L., Simmonds, T.J. et al. HiFiAdapterFilt, a memory efficient read processing pipeline, prevents occurrence of adapter sequence in PacBio HiFi reads and their negative impacts on genome assembly. BMC Genomics 23, 157 (2022). https://doi.org/10.1186/s12864-022-08375-1"; exit 0 ;;
v | version ) echo $version; exit 0 ;;
h | help ) echo $USAGE; exit 0;;
?*) echo "Illegal option -$OPT"; echo $USAGE; exit 2;;
esac
done

if ((OPTIND == 1))
then
    printf "No options specified. \nFiltering files in working directory. \n"
fi

shift $((OPTIND - 1)) # remove parsed options and args from $@ list

## Set variables based on options

if [ -f ${prefix}.temp_file_list ]
then
	rm ${prefix}.temp_file_list
fi

if [ ${prefix} = "all" ]
then
        ls | grep "bam\|fq\|fastq" >${prefix}.temp_file_list

else
	if ls ${prefix}*.f*q* >/dev/null 2>&1
	then
		for x in `ls ${prefix}*.f*q*`; do echo $x >>${prefix}.temp_file_list; done
	fi

	if ls ${prefix}*.bam >/dev/null 2>&1
	then
		for x in `ls ${prefix}*.bam`; do echo $x >>${prefix}.temp_file_list; done
	fi
fi


##

reads_pref=$(for x in `cat ${prefix}.temp_file_list`; do ls $x | sed 's/\.fastq.gz//' | sed 's/\.fq.gz//' | sed 's/\.bam//' | sed 's/\.fastq//' | sed 's/\.fq//'; done)
read_path=$(dirname ${prefix}*)
read_path_str=$(echo ${read_path} | cut -d" " -f 1)

## Create out directory if necessary

if [ ! -d ${outdir} ]
then 
	mkdir ${outdir}
fi

for x in `cat ${prefix}.temp_file_list`; do echo $x; done

if [ -z ${reads_pref} ]
then
	echo "No files of proper name or filetype recognized. Exiting $(date)."
	exit
fi

## BLAST raw files

for x in `echo ${reads_pref}`
do
	touch ${outdir}/${x}.stats
	echo "Started on $(date)" >>${outdir}/${x}.stats

	if [ -s ${read_path_str}/${x}.fastq.gz ]
	then

	   if command -v pigz &> /dev/null
	   then
		pigz -cd -p ${threads} ${read_path_str}/${x}.fastq.gz > ${read_path_str}/${x}.fastq 
	   else
	   	zcat ${read_path_str}/${x}.fastq.gz > ${read_path_str}/${x}.fastq 
	   fi

	   echo "Identifying reads in ${x}.fastq.gz with adapter contamination on $(date)."
	   cat ${read_path_str}/${x}.fastq | sed -n '1~4s/^@/>/p;2~4p' > ${read_path_str}/${x}.fasta
	   blastn -db $DBpath/pacbio_vectors_db -query ${read_path_str}/${x}.fasta -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 > ${outdir}/${x}.contaminant.blastout &
	   wait
	   echo "Creating blocklist for ${x}.fastq.gz on $(date)."
	   cat ${outdir}/${x}.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' -v var1="${adapterlength}" -v var2="${pctmatch}" '{if (($2 ~ /NGB00972/ && $3 >= var2 && $4 >= var1) || ($2 ~ /NGB00973/ && $3 >= 97 && $4 >= 34)) print $1}' | sort -u > ${outdir}/${x}.blocklist &
	   wait
	   echo "Removing adapter contaminated reads from ${x}.fastq.gz on $(date)." 

	   if command -v pigz &> /dev/null
	   then
	   	cat ${read_path_str}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" | pigz -p ${threads} --fast > ${outdir}/${x}.filt.fastq.gz
	   else
	   	cat ${read_path_str}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" | gzip -1 > ${outdir}/${x}.filt.fastq.gz
	   fi

	   f=`cat ${outdir}/${x}.blocklist | wc -l` #number of adapter contaminated reads
	   r1=`zcat ${read_path_str}/${x}.fastq.gz | wc -l` 
	   r2=`awk -v r1=$r1 'BEGIN{ans=r1/4; print ans}'` #number of ccs reads
	   p1=`awk -v n1=$r2 -v n2=$f 'BEGIN{ans=n2/n1*100; print ans}'` #proportion of adapter contaminated reads
	   r3=`awk -v r2=$r2 -v f=$f 'BEGIN{ans=r2-f; print ans}'` #number of reads retained
	   p2=`awk -v p1=$p1 'BEGIN{ans=100-p1; print ans}'` #proportion of reads retained
	   echo "For the" ${x} "dataset:" >>${outdir}/${x}.stats 
	   echo "Removing reads containing adapters a minimum of ${adapterlength} bp in length and ${pctmatch}% match." >>${outdir}/${x}.stats
	   echo "" >>${outdir}/${x}.stats
	   echo "Number of ccs reads:" $r2 >>${outdir}/${x}.stats
	   echo "Number of adapter contaminated ccs reads:" $f "("$p1"% of total)" >>${outdir}/${x}.stats
	   echo "Number of ccs reads retained:" $r3 "("$p2"% of total)" >>${outdir}/${x}.stats
	   echo "" >>${outdir}/${x}.stats
	   echo "Finished on $(date)" >>${outdir}/${x}.stats
	   rm ${read_path_str}/${x}.fastq
	   rm ${read_path_str}/${x}.fasta

	elif [ -s ${read_path_str}/${x}.fq.gz ]
	then
	   echo "Identifying reads in ${x}.fq.gz with adapter contamination on $(date)."

	   if command -v pigz &> /dev/null
	   then
		pigz -cd -p ${threads} ${read_path_str}/${x}.fq.gz > ${read_path_str}/${x}.fq
	   else
	   	zcat ${read_path_str}/${x}.fq.gz > ${read_path_str}/${x}.fq 
	   fi

	   cat ${read_path_str}/${x}.fq | sed -n '1~4s/^@/>/p;2~4p' > ${read_path_str}/${x}.fasta 
	   blastn -db $DBpath/pacbio_vectors_db -query ${read_path_str}/${x}.fasta -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 > ${outdir}/${x}.contaminant.blastout &
	   wait
	   echo "Creating blocklist for ${x}.fq.gz on $(date)."
	   cat ${outdir}/${x}.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' -v var1="${adapterlength}" -v var2="${pctmatch}" '{if (($2 ~ /NGB00972/ && $3 >= var2 && $4 >= var1) || ($2 ~ /NGB00973/ && $3 >= 97 && $4 >= 34)) print $1}' | sort -u > ${outdir}/${x}.blocklist &
	   wait
	   echo "Removing adapter contaminated reads from ${x}.fastq.gz on $(date)." 

	   if command -v pigz &> /dev/null
	   then
	   	cat ${read_path_str}/${x}.fq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" | pigz -p ${threads} --fast > ${outdir}/${x}.filt.fastq.gz
	   else
	   	cat ${read_path_str}/${x}.fq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" | gzip -1 > ${outdir}/${x}.filt.fastq.gz
	   fi

	   f=`cat ${outdir}/${x}.blocklist | wc -l` #number of adapter contaminated reads
	   r1=`zcat ${read_path_str}/${x}.fq.gz | wc -l` 
	   r2=`awk -v r1=$r1 'BEGIN{ans=r1/4; print ans}'` #number of ccs reads
	   p1=`awk -v n1=$r2 -v n2=$f 'BEGIN{ans=n2/n1*100; print ans}'` #proportion of adapter contaminated reads
	   r3=`awk -v r2=$r2 -v f=$f 'BEGIN{ans=r2-f; print ans}'` #number of reads retained
	   p2=`awk -v p1=$p1 'BEGIN{ans=100-p1; print ans}'` #proportion of reads retained
	   echo "For the" ${x} "dataset:" >>${outdir}/${x}.stats 
	   echo "Removing reads containing adapters a minimum of ${adapterlength} bp in length and ${pctmatch}% match." >>${outdir}/${x}.stats
	   echo "" >>${outdir}/${x}.stats
	   echo "Number of ccs reads:" $r2 >>${outdir}/${x}.stats
	   echo "Number of adapter contaminated ccs reads:" $f "("$p1"% of total)" >>${outdir}/${x}.stats
	   echo "Number of ccs reads retained:" $r3 "("$p2"% of total)" >>${outdir}/${x}.stats
	   echo "" >>${outdir}/${x}.stats
	   echo "Finished on $(date)" >>${outdir}/${x}.stats
	   rm ${read_path_str}/${x}.fq
	   rm ${read_path_str}/${x}.fasta
	elif [ -s ${read_path_str}/${x}.bam ]
	then

	   if [ ! -s ${read_path_str}/${x}.fastq ]
		then
    		echo "Converting ${x}.bam to ${x}.fastq on $(date)"
    		bamtools convert -format fastq -in ${read_path_str}/${x}.bam -out ${read_path_str}/${x}.fastq &
	   fi

	   if [ ! -s ${read_path_str}/${x}.fasta ]
	   then
        	echo "Converting ${x}.bam to ${x}.fasta on $(date)."    
        	bamtools convert -format fasta -in ${read_path_str}/${x}.bam -out ${read_path_str}/${x}.fasta &
        	wait
	   fi

	   echo "Identifying reads in ${x}.bam with adapter contamination on $(date)."
	   blastn -db $DBpath/pacbio_vectors_db -query ${outdir}/${x}.fasta -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 > ${outdir}/${x}.contaminant.blastout &
	   wait
	   echo "Creating blocklist for ${x}.bam on $(date)."
	   cat ${outdir}/${x}.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' -v var1="${adapterlength}" -v var2="${pctmatch}" '{if (($2 ~ /NGB00972/ && $3 >= var2 && $4 >= var1) || ($2 ~ /NGB00973/ && $3 >= 97 && $4 >= 34)) print $1}' | sort -u > ${outdir}/${x}.blocklist &
	   wait
	   echo "Removing adapter contaminated reads from ${x}.bam on $(date)." 

	   if command -v pigz &> /dev/null
	   then
	   	cat ${read_path_str}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" | pigz -p ${threads} --fast > ${outdir}/${x}.filt.fastq.gz
	   else
	   	cat ${read_path_str}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" | gzip -1 > ${outdir}/${x}.filt.fastq.gz
	   fi

	   f=`cat ${outdir}/${x}.blocklist | wc -l` #number of adapter contaminated reads
	   r1=`cat ${outdir}/${x}.fastq | wc -l` 
	   r2=`awk -v r1=$r1 'BEGIN{ans=r1/4; print ans}'` #number of ccs reads
	   p1=`awk -v n1=$r2 -v n2=$f 'BEGIN{ans=n2/n1*100; print ans}'` #proportion of adapter contaminated reads
	   r3=`awk -v r2=$r2 -v f=$f 'BEGIN{ans=r2-f; print ans}'` #number of reads retained
	   p2=`awk -v p1=$p1 'BEGIN{ans=100-p1; print ans}'` #proportion of reads retained
	   echo "For the" ${x} "dataset:" >>${outdir}/${x}.stats 
	   echo "Removing reads containing adapters a minimum of ${adapterlength} bp in length and ${pctmatch}% match." >>${outdir}/${x}.stats
	   echo "" >>${outdir}/${x}.stats
	   echo "Number of ccs reads:" $r2 >>${outdir}/${x}.stats
	   echo "Number of adapter contaminated ccs reads:" $f "("$p1"% of total)" >>${outdir}/${x}.stats
	   echo "Number of ccs reads retained:" $r3 "("$p2"% of total)" >>${outdir}/${x}.stats
	   echo "" >>${outdir}/${x}.stats
	   echo "Finished on $(date)" >>${outdir}/${x}.stats
	   rm ${read_path_str}/${x}.fasta
	   rm ${read_path_str}/${x}.fastq

	elif [ -s ${read_path_str}/${x}.fastq ]
	then
	   echo "Identifying reads in ${x}.fastq with adapter contamination on $(date)."
	   cat ${read_path_str}/${x}.fastq | sed -n '1~4s/^@/>/p;2~4p' > ${read_path_str}/${x}.fasta
	   blastn -db $DBpath/pacbio_vectors_db -query ${read_path_str}/${x}.fasta -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 > ${outdir}/${x}.contaminant.blastout &
	   wait
	   echo "Creating blocklist for ${x}.fastq on $(date)."
	   cat ${outdir}/${x}.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' -v var1="${adapterlength}" -v var2="${pctmatch}" '{if (($2 ~ /NGB00972/ && $3 >= var2 && $4 >= var1) || ($2 ~ /NGB00973/ && $3 >= 97 && $4 >= 34)) print $1}' | sort -u > ${outdir}/${x}.blocklist &
	   wait
	   echo "Removing adapter contaminated reads from ${x}.fastq on $(date)."

	   if command -v pigz &> /dev/null
	   then
	   	cat ${read_path_str}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" | pigz -p ${threads} --fast > ${outdir}/${x}.filt.fastq.gz
	   else
	   	cat ${read_path_str}/${x}.fastq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" | gzip -1 > ${outdir}/${x}.filt.fastq.gz
	   fi

	   f=`cat ${outdir}/${x}.blocklist | wc -l` #number of adapter contaminated reads
	   r1=`cat ${read_path_str}/${x}.fastq | wc -l` 
	   r2=`awk -v r1=$r1 'BEGIN{ans=r1/4; print ans}'` #number of ccs reads
	   p1=`awk -v n1=$r2 -v n2=$f 'BEGIN{ans=n2/n1*100; print ans}'` #proportion of adapter contaminated reads
	   r3=`awk -v r2=$r2 -v f=$f 'BEGIN{ans=r2-f; print ans}'` #number of reads retained
	   p2=`awk -v p1=$p1 'BEGIN{ans=100-p1; print ans}'` #proportion of reads retained
	   echo "For the" ${x} "dataset:" >>${outdir}/${x}.stats 
	   echo "Removing reads containing adapters a minimum of ${adapterlength} bp in length and ${pctmatch}% match." >>${outdir}/${x}.stats
	   echo "" >>${outdir}/${x}.stats
	   echo "Number of ccs reads:" $r2 >>${outdir}/${x}.stats
	   echo "Number of adapter contaminated ccs reads:" $f "("$p1"% of total)" >>${outdir}/${x}.stats
	   echo "Number of ccs reads retained:" $r3 "("$p2"% of total)" >>${outdir}/${x}.stats
	   echo "" >>${outdir}/${x}.stats
	   echo "Finished on $(date)" >>${outdir}/${x}.stats
	   rm ${read_path_str}/${x}.fasta
	elif [ -s ${read_path_str}/${x}.fq ]
	then
	   echo "Identifying reads in ${x}.fq with adapter contamination on $(date)."
	   cat ${read_path_str}/${x}.fq | sed -n '1~4s/^@/>/p;2~4p' > ${read_path_str}/${x}.fasta
	   blastn -db $DBpath/pacbio_vectors_db -query ${read_path_str}/${x}.fasta -num_threads ${threads} -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust no -soft_masking true -evalue 700 -searchsp 1750000000000 -outfmt 6 > ${outdir}/${x}.contaminant.blastout &
	   wait
	   echo "Creating blocklist for ${x}.fq on $(date)."
	   cat ${outdir}/${x}.contaminant.blastout | grep 'NGB0097' | awk -v OFS='\t' -v var1="${adapterlength}" -v var2="${pctmatch}" '{if (($2 ~ /NGB00972/ && $3 >= var2 && $4 >= var1) || ($2 ~ /NGB00973/ && $3 >= 97 && $4 >= 34)) print $1}' | sort -u > ${outdir}/${x}.blocklist &
	   wait
	   echo "Removing adapter contaminated reads from ${x}.fq on $(date)."

	   if command -v pigz &> /dev/null
	   then
	   	cat ${read_path_str}/${x}.fq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" | pigz -p ${threads} --fast > ${outdir}/${x}.filt.fastq.gz
	   else
	   	cat ${read_path_str}/${x}.fq | paste - - - - | grep -v -f ${outdir}/${x}.blocklist -F | tr "\t" "\n" | gzip -1 > ${outdir}/${x}.filt.fastq.gz
	   fi

	   f=`cat ${outdir}/${x}.blocklist | wc -l` #number of adapter contaminated reads
	   r1=`cat ${read_path_str}/${x}.fq | wc -l` 
	   r2=`awk -v r1=$r1 'BEGIN{ans=r1/4; print ans}'` #number of ccs reads
	   p1=`awk -v n1=$r2 -v n2=$f 'BEGIN{ans=n2/n1*100; print ans}'` #proportion of adapter contaminated reads
	   r3=`awk -v r2=$r2 -v f=$f 'BEGIN{ans=r2-f; print ans}'` #number of reads retained
	   p2=`awk -v p1=$p1 'BEGIN{ans=100-p1; print ans}'` #proportion of reads retained
	   echo "For the" ${x} "dataset:" >>${outdir}/${x}.stats 
	   echo "Removing reads containing adapters a minimum of ${adapterlength} bp in length and ${pctmatch}% match." >>${outdir}/${x}.stats
	   echo "" >>${outdir}/${x}.stats
	   echo "Number of ccs reads:" $r2 >>${outdir}/${x}.stats
	   echo "Number of adapter contaminated ccs reads:" $f "("$p1"% of total)" >>${outdir}/${x}.stats
	   echo "Number of ccs reads retained:" $r3 "("$p2"% of total)" >>${outdir}/${x}.stats
	   echo "" >>${outdir}/${x}.stats
	   echo "Finished on $(date)" >>${outdir}/${x}.stats
	   rm ${read_path_str}/${x}.fasta
	fi
done

if [ -f ${prefix}.temp_file_list ]
then
	rm ${prefix}.temp_file_list
fi
