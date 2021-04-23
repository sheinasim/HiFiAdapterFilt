# HiFiAdapterFilt
Convert .bam to .fastq and remove reads with remnant PacBio adapter sequences  

Dependencies:

* Bamtools 
* BLAST+

Add script and database to your path using:  

> export PATH=$PATH:[PATH TO HiFiAdapterFilt]  
> export PATH=$PATH:[PATH TO HiFiAdapterFilt]/DB  

Usage:  
> sh pbadapterfilt.sh [ -b .bam prefix ] [ -t Number of threads for blastn. Default=1] [ -o outdirectory prefix Default=. ]  

If no arguments are provided, the script will run on all .bam files in the working directory.

### If this script is useful to you, please cite the following in your publication:

@software{HiFiAdapterFilt,
  author = {{Sheina B. Sim}},
  title = {HiFiAdapterFilt},
  url = {https://hadoop.apache.org},
  DOI = {10.5281/zenodo.4716418},
  version = {v1.0.0},
  date = {2021-04-23},
}

<a href="https://doi.org/10.5281/zenodo.4716418"><img src="https://github.com/sheinasim/HiFiAdapterFilt/blob/master/zenodo.4716418.png" width="250" title="HiFiAdapterFilt DOI" alt="DOI"></a>

Sheina B. Sim  
USDA-ARS  
US Pacific Basin Agricultural Research Service  
Hilo, Hawaii, 96720 USA  
sheina.sim@usda.gov  

This script is in the public domain in the United States per 17 U.S.C. § 105
