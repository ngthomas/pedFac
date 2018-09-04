# pedfac
A factor graph based Markov chain Monte Carlo (MCMC) pedigree sampler

pedFac.py is a wrapper script that oversees the complete workflow of the MCMC based pedigree sampler.

## Requirements and Installation

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

This Python script expects the user to provide the path to a valid space-separated genotype file named "genotype.txt" so that it can generate all necessary secondary files for the C-script pedigree sampler to run. Once the sampling iterations are completed, it returns a summary of the MCMC output and also details of the sampled pedigrees.  

### The genotype file:

The genotype file is a space-separated file that contains genotype and metadata information
from the individuals whose pedigree is to be determined.   
Each line in the file holds information about a single individual ordered as follows: 

- unique individual id | is the indiv observed? | sex of individual | birth year | genotype(s) information
  
For example, for an observed male individual who

- is born in the middle of 1986  
- with heterozygous allele of "1" and "0" for the first locus and  
- homogygous allele of "3" for the second locus,  

the entry would look like this: `10 1 1 1986.5 1 0 3 3`   

A more detailed breakdown for each column/field in each individual row in the file 
is given here:

* **column 1**: unique individual identifer, it must be a positive, unique integer for each row  
* **column 2**: an integer flag saying whether the individual is observed (i.e. whether genotype
data are available for it)  Must be 1 or 0,  corresponding to yes or no.  
* **column 3**: sex of the individual: 0, 1, or 2. - corresponding to unknown, male or female.
* **column 4**: birth year: a real number value allowing a decimal point,  e.g. 1994.10  
* **column 5+6 (7+8) ... (x+(x+1))**: columns holding the genotype information of a 'diploid organism'. Each locus occupies
two columns (one for each gene copy). The allele information must be giving in one of the following forms:  
    - *As an integer*. If genotype information is not known, use -1. e.g -1 -1. If data are from biallelic SNPs, we recommend using 0
    for one allele and 1 for the other.  
    - *string for both alleles for any number of loci*. For string based genotype i.e. haplotype, this program recognizes standard A, T, C, and G bases, e.g ACAAT ATCAA. If the genotype info is not known, use N e.g. `N N`. (NOT CLEAR THOMAS, IS THIS REALLY FOR HAPLOTYPES? DO WE USE JUST ONE COLUMN FOR ALL MARKERS? NEED TO TALK THIS OVER).  
    - *comma-separated genotype classes and associated posterior probabilties*. (TECHNICALLY THESE ARE NORMALIZED LIKELIHOODS, NOT POSTERIOR PROBS, NO?) For now, we only accept the case of biallelic genotypes with 3 possible genotypic classes - 0, 1, 2.  `0` represents a homozygote for the major allele, `2` represents a homozygote for the minor allele, and `1` represents a heterozygote.  There is one column for each locus, and it is a string of comma-separate genotype classes (i.e. 0,1,2) followed by a string of their respective genotype probabilities.  For example `0,1,2 0.9,0.3,0.2`   
    - (This feature is not available yet --) For dealing with multiallelic classes (not avail yet), you will be only need to provide the 
    top four (if any) genotype classes. Any remaining prob will be spread uniformly over the other unlisted
    classes, e.g. `AA|AT,AC|AT,AA|AC,AT|AT 0.5,0.1,0.1,0.1`
  
### Optional marker input file: "marker_info.txt":  

This space-separated input file -- "marker_info.txt" --- can only be given for biallelic marker data.  It holds
metadata information regarding the status of the genotype markers.  

The first line of the file begins with the tag "name" followed by a space-separated list of the names of the loci.

In subsequent lines, one can choose to put genotyping error and allele frequency information by starting the lines with
the appropriate tags, as listed below:

- `gerror` :  genotyping error rate at the locus. A positive value from 0 to 1. It describes rates of genotyping error from a very 
simplistic model, the  "alpha" model: with probability of (1-a/n) the true genotype is the observed genotype and
with probability a/n the genotype is drawn randomly from the population
genotype frequencies. If genotype posterior values are given, the `gerror` information is ignored.   
- `afreq`: allele frequency of the "1" allele (as opposed to the "0") allele. Must be a positive float value from 0 to 1.   

e.g., the file could look like this:
```
name SNP_1 SNP_2 SNP_3 SNP_4  
gerror 0.2 0.2 0.2 0.2  
afreq 0.01 0.04 0.1 0.23  
```

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

