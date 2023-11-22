# HiFiAdapterFilt
Convert .bam to .fastq and remove reads with remnant PacBio adapter sequences  

Dependencies:

* BamTools 
* BLAST+

Optional:

* NCBI FCS Adaptor
* pigz

Add script and database to your path using:  

> ``export PATH=$PATH:[PATH TO HiFiAdapterFilt]``  
> ``export PATH=$PATH:[PATH TO HiFiAdapterFilt]/DB``  

### Usage
  
> ``bash hifiadapterfilt.sh [ -p file Prefix ] [ -l minimum Length of adapter match to remove. Default=44 ] [ -m minimum percent Match of adapter to remove. Default=97 ] [ -t Number of threads for blastn. Default=8 ] [ -o outdirectory prefix Default=. ]`` 

All flags are optional. 

If no -p argument is provided, the script will run on all sequence files (.bam, .fastq, .fastq.gz, .fq, .fq.gz) in the working directory.

If using FCS adaptor to detect adapter contaminated reads use the hifiadapterfiltFCS.sh script

> ``bash hifiadapterfiltFCS.sh -f <FCS adaptor output file> -r <HiFi reads file> [-t Number of threads for pigz. Defualt=8] [-o outdirectory prefix Default=.]``

## Outputs

* {prefix}.contaminant.blastout (Output of BLAST search)
* {prefix}.blocklist (Headers of PB adapter contaminated reads to be removed)
* {prefix}.filt.fastq.gz (Fastq reads free of PB adapter sequence ready for assembly)
* {prefix}.stats (File with simple math on number of reads removed, etc)

### Citation

If this script is useful to you, please cite the following in your publication:

```
@article{HiFiAdapterFilt,
   author = {Sim, Sheina B. and Corpuz, Renee L. and Simmonds, Tyler J. and Geib, Scott M.},
   title = {HiFiAdapterFilt, a memory efficient read processing pipeline, prevents occurrence of adapter sequence in PacBio HiFi reads and their negative impacts on genome assembly},
   journal = {BMC Genomics},
   volume = {23},
   number = {1},
   pages = {157},
   ISSN = {1471-2164},
   DOI = {10.1186/s12864-022-08375-1},
   url = {https://doi.org/10.1186/s12864-022-08375-1},
   year = {2022},
   type = {Journal Article}
}

```

Sheina B. Sim  
USDA-ARS  
US Pacific Basin Agricultural Research Service  
Hilo, Hawaii, 96720 USA  
sheina.sim@usda.gov  

This script is in the public domain in the United States per 17 U.S.C. § 105
