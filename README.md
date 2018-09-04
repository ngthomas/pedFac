# pedfac
A factor graph based Markov chain Monte Carlo (MCMC) pedigree sampler

pedFac.py is a wrapper script that oversees the complete workflow of the MCMC based pedigree sampler.

## Requirement

This package requires python v 3.0+ to be installed in either a Linux or a Mac OS environment. Check out : [conda.io/miniconda.html](conda.io/miniconda.html) to install or update your local python.    

After python v 3.0+ is installed, the next step is to ensure that the python package `numpy` is installed.  Do so with

```sh
conda install numpy
```

Once python and its requisite packages are installed, you can install the pedfac software by cloning this
GitHub repoistory, like so:

```sh
git clone https://github.com/ngthomas/pedfac.git
```

You can perform a test run as follows: 

```sh
cd pedfac
python bin/run-pedfac -i example/case_0/ -n 5
```

## About run-pedfac

This Python script expects the user to provide the path to a valid space separate genotype file named "genotype.txt" so that it can generate all necessary secondary files for the C-script pedigree sampler to run. Once the sampling iterations are completed, it returns the summary and also the details of the sampled pedigrees.  

### About the genotype file:
The genotype file is a space separate file that contains individual's genotype and meta information.   
Each row is an individual entry with its associate genotype information, in the order as follows: unique indiv ID | is the indiv observed? | sex of individual | birth year | genotype(s) information.    
  
For example, for an observed male indiv who   
    - is born in the middle of 1986  
    - with heterozygous allele of "1" and "0" for a locus and  
    - homogygous allele of "3" for another locus,  
the entry would look like this: 10 1 1 1986.5 1 0 3 3   

This is more of a detailed breakdown for each column field:  
column 1: unique indiv ID - a positive unique integer for each row  
column 2: is the individual observed? 1 or 0 - corresponds to yes or no.  
column 3: sex of the individual: 0, 1, 2. - corresponds to unknown, male or female   
column 4: birth year: flow value. allow decimal value e.g. 1994.10  
column 5+6(7+8)....(x+(x+1)): genotype information of a 'diploid organism'. Must provide the alleles info in one of the following forms:  
    - integer. If genotype information is not known, use -1. e.g -1 -1. If it is a biallelic SNP, we recommend using 0,1.  
    - string for both alleles for any number of loci. For string based genotype i.e. haplotype, this program recognizes standard A, T, C, and G bases, e.g ACAAT ATCAA. If the genotype info is not known, use N e.g. N N.  
    - comma-separate genotype class and posterior probabilty. For now, we only accept biallelic genotype case with 3 possible classes - 0, 1, 2. Let 0 or 2 be the homozygous case of possessing the more common or more rare allele, respectively. Let 1 be the heterozygous case. The first column is a string of comma-separate genotype classes (i.e. 0,1,2) followed by a column of their respective genotype probability.  
    e.g 0,1,2 0.9,0.3,0.2   
    (This feature is not available yet --) As for the case of dealing with multiallelic classes (not avail yet), you will be only need to provide the top four (if any) genotype classes. Any remaining prob will be spread uniformly for other unlisted classes,   
    e.g. AA|AT,AC|AT,AA|AC,AT|AT 0.5,0.1,0.1,0.1   
  
Optional marker input file (only for biallelic markers) - "marker_info.txt":  
This space separate input file -- "marker_info.txt" --- holds meta information regarding the status of the genotype markers.  

The first line of the file begins with the tag "name" following by names of locus. (space separated)  
One can choose any of the following description (if any) associated with the markers:  

- gerror :  genotype error. positive value from 0 to 1. It describes the rates the genotype error with the most simplistic model - alpha model. prob of (1-a/n) for observed genotype class. If genotype posterior value is already reported, this info is going to be ignored.   

- afreq: allelic frequency of the alternative allele - "1", as opposed to "0". Positive float value from 0 to 1.   
  
e.g.
name SNP_1 SNP_2 SNP_3 SNP_4  
gerror 0.2 0.2 0.2 0.2  
afreq 0.01 0.04 0.1 0.23  
  
### About this wrapper parameters:

    -i/--inputPath: directory path that contains the input "geno.txt" and optional "marker_info.txt". String. Required.
    -o/--outputPath: directory path to store intermed and final output files. String. Optional. If not specified, use input path

    Regarding sampling:
    -r/--randomSeed: random seed to pass on sampler. Positive integer; a randomly generated value as default
    -n/--nIter: number of sampling iteration. Positive integer; 1 as default
    -c/--cyclicChoice: choices of handling loops. [0 (default), 1, or 2].
        0 - not allowing loops;
        1 - throttle method;
        2 - decimation method

    -f/--observeFrac: assumed sampling fraction. Float value from 0 to 1; 0.8 as default. However, if you don't want to impose any prior knowledge about sampling fraction, use the value of -1.

    -u/--maxUnobs: maximum number of unobserved individuals allowed in between any two individuals. Nonnegative integer; 1 as default
    -m/--maxGen: number of predecessor generation(s) considered beyond the earliest observed generation. Nonnegative integer; 0 as default. Setting it as 0 means that individuals of the earliest observed generation are treated as founders.

    Regarding specie life history:
    -s/--minAge: minimal age of sexual maturation or fecundity (in year). Positive float value; 1 as default
    -a/--maxAge: maximum age of sexual maturation or fecundity (in year). Positive float value; 1 as default

    Regarding genotype marker: (need to generate intermediate summary files of how markers are compressed/collapsed )
    -hm/--haploMethod: Selected method in the case of handling multiallelic markers. Positive integer 0 - 2; 0 as default.
        0 - taking the most informative allele whose frequency is closest to 0.5
        1 - (not avail) deconstructing haplotype into a set of nucleotide units
        2 - (not avail) reduce the multiallelic basis into n class of binomial switches
    -g/--genoerr: Assumed background genotype error rate in the form of epsilon. Float value from 0 to 1; 0.02 as default. If the genotype error row - 'gerror' of marker_info.txt is provided, this param will be overridden.

