# HiFiAdapterFilt
Convert .bam to .fastq and remove reads with remnant PacBio adapter sequences  

Dependencies:

* Bamtools 
* BLAST+

Add script and database to your path using:  

> export PATH=$PATH:[PATH TO HiFiAdapterFilt]  
> export PATH=$PATH:[PATH TO HiFiAdapterFilt]/DB  

### Usage
  
> bash pbadapterfilt.sh [ -p file Prefix ] [ -l minimum Length of adapter match to remove. Default=44 ] [ -m minimum percent Match of adapter to remove. Default=97 ] [ -t Number of threads for blastn. Default=8 ] [ -o outdirectory prefix Default=. ]  

All flags are optional. 

If no -p argument is provided, the script will run on all sequence files (.bam, .fastq, .fastq.gz, .fq, .fq.gz) in the working directory.

## Outputs

* {prefix}.contaminant.blastout (Output of BLAST search)
* {prefix}.blocklist (Headers of PB adapter contaminated reads to be removed)
* {prefix}.filt.fastq.gz (Fastq reads free of PB adapter sequence ready for assembly)
* {prefix}.stats (File with simple math on number of reads removed, etc)

### Citation

If this script is useful to you, please cite the following in your publication:

```
@software{HiFiAdapterFilt,
  author = {Sim, Sheina B.},
  title = {HiFiAdapterFilt},
  url = {https://github.com/sheinasim/HiFiAdapterFilt},
  DOI = {10.5281/zenodo.4716418},
  version = {v1.0.0},
  date = {2021-04-23},
}
```

<a href="https://doi.org/10.5281/zenodo.4716418"><img src="https://github.com/sheinasim/HiFiAdapterFilt/blob/master/zenodo.4716418.png" width="250" title="HiFiAdapterFilt DOI" alt="DOI"></a>

Sheina B. Sim  
USDA-ARS  
US Pacific Basin Agricultural Research Service  
Hilo, Hawaii, 96720 USA  
sheina.sim@usda.gov  

This script is in the public domain in the United States per 17 U.S.C. § 105
