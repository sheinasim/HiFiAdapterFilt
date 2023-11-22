#!/bin/bash

version=3.0.1
threads=8
outdir=$(pwd)

# Function to display script usage
usage() {
    echo "Usage: ./hifiadapterfiltFCS.sh -f <FCS adapter output file (fcs_adaptor_report.txt)> -r <HiFi reads file> [-t] [-o] [-c] [-v]"
    echo "Options:"
    echo "  -f <fcs file>: Required FCS adapter output file (fcs_adaptor_report.txt)"
    echo "  -r <sequence file>: Required equence file"
    echo "  -t: Optional number of threads (Default=8)"
    echo "  -o: Optional output directory (Default=.)"
    echo "  -c: Optional print citation"
    echo "  -v: Optional print version"
    exit 1
}

# Check if no arguments are passed
if [ $# -eq 0 ]; then
    usage
fi

# Parsing the arguments using a while loop
while getopts ":f:r:t:o:cv" opt; do
    case $opt in
        f)
            fcs=$OPTARG
            ;;
        r)
            reads=$OPTARG
            ;;
        t)
            threads=$OPTARG
            ;;
        o)
            outdir=$OPTARG
            ;;
        c)
            echo "Sim, S.B., Corpuz, R.L., Simmonds, T.J. et al. HiFiAdapterFilt, a memory efficient read processing pipeline, prevents occurrence of adapter sequence in PacBio HiFi reads and their negative impacts on genome assembly. BMC Genomics 23, 157 (2022). https://doi.org/10.1186/s12864-022-08375-1"
            exit 0 
            ;;
        v)
            echo $version
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            usage
            ;;
    esac
done

# Check if the required FCS adaptor output file is provided
if [ -z "$fcs" ] || [ -z "$reads" ]; then
    echo "A file must be provided using -f option."
    usage
fi

# Displaying the values after parsing the arguments
echo "FCS output File: $fcs"
echo "Reads: $reads"
echo "Threads: $threads"
echo "Outdir: $outdir"

## Create out directory if necessary

if [ ! -d "${outdir}" ]
then
        mkdir ${outdir}
fi

if [ -z "${reads}" ] || [ -z "${fcs}" ]
then
        echo "One or both required files are missing. Exiting $(date)."
        exit
fi

filetype=`ls ${reads} | awk -F "." '{if ($NF == "gz") print $(NF-1)"."$NF; else print $NF}'`
readsprefix=`ls $reads | sed "s/\.${filetype}//"`

if [ $filetype == "fq.gz" ]
then
    mv ${readsprefix}.fq.gz ${readsprefix}.fastq.gz
fi

if [ $filetype == "bam" ]
then
    bamtools convert -format fastq -in ${readsprefix}.bam -out ${readsprefix}.fastq
    bamtools convert -format fasta -in ${readsprefix}.bam -out ${readsprefix}.fasta
elif [ $filetype == "fastq.gz" ]
then
    if command -v pigz &> /dev/null
    then
        pigz -cd -p ${threads} ${readsprefix}.fastq.gz > ${readsprefix}.fastq
    else
        zcat ${readsprefix}.fastq.gz > ${readsprefix}.fastq
    fi
fi

cat ${fcs} | grep "NGB00972.1" | awk -v OFS='\t' '{print $1}' | sort -u > ${readsprefix}.blocklist

if command -v pigz &> /dev/null
then
    cat ${readsprefix}.fastq | paste - - - - | grep -v -f ${readsprefix}.blocklist -F | tr "\t" "\n" | pigz -p 40 --fast > ${readsprefix}.fcsfilt.fastq.gz
else
    cat ${readsprefix}.fastq | paste - - - - | grep -v -f ${readsprefix}.blocklist -F | tr "\t" "\n" | gzip - >${readsprefix}.fcsfilt.fastq.gz
fi

f=`cat ${readssprefix}.blocklist | wc -l` #number of adapter contaminated 
r1=`cat ${readssprefix}.fastq | wc -l`
r2=`awk -v r1=$r1 'BEGIN{ans=r1/4; print ans}'` #number of ccs reads
p1=`awk -v n1=$r2 -v n2=$f 'BEGIN{ans=n2/n1*100; print ans}'` #proportion of adapter contaminated reads
r3=`awk -v r2=$r2 -v f=$f 'BEGIN{ans=r2-f; print ans}'` #number of reads retained
p2=`awk -v p1=$p1 'BEGIN{ans=100-p1; print ans}'` #proportion of reads retained
echo "For the" ${readssprefix} "dataset:" >>${readssprefix}.stats
echo "" >>${readssprefix}.stats
echo "Number of ccs reads:" $r2 >>${readssprefix}.stats
echo "Number of adapter contaminated ccs reads:" $f "("$p1"% of total)" >>${readssprefix}.stats
echo "Number of ccs reads retained:" $r3 "("$p2"% of total)" >>${readssprefix}.stats
echo "" >>${readssprefix}.stats
echo "Finished on $(date)" >>${readssprefix}.stats
done
