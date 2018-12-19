# BetaScan
This software implements the Beta statistic to detect ancient balancing selection. Our paper describing this statistic and its application is available [here](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx209/3988103/Detecting-Long-term-Balancing-Selection-using). BetaScan takes in a file of variant positions and frequencies, and calculates Beta in a sliding window approach. It allows the user to choose appropriate parameter choices, and outputs the value of Beta for each variant.

Any feedback or questions are very welcome. You can e-mail Katie at ksiewert@pennmedicine.upenn.edu or post a github issue. We know that programs written by other people can be difficult to use, so we’ve tried our best to make this program simple and intuitive. That being said, bioinformatics is bioinformatics, and issues will arise, so don’t hesitate to contact us!

# 1000 Genomes Beta Scores
If you would like the Beta Scores for each population in the 1000 Genomes dataset, they are available [here](http://coruscant.itmat.upenn.edu/data/SiewertEA_Full_BetaScores.tar.gz) (warning: this is a 1.8 GB gzipped file). If you just want to look at the top scoring haplotypes in each population, that data is available [here](http://coruscant.itmat.upenn.edu/data/SiewertEA_BetaScores.tar.gz).

# Recent Updates
12/11/18: BetaScan version 2 is released! This version has two main improvements: it implements a new Beta statistic which uses substitutions, and it can now standardize each Beta statistic by its variance. These updates are detailed in our [bioRxiv preprint](https://www.biorxiv.org/content/early/2018/12/17/497255).

10/15/18: You can now specify the output file name using the -o flag. Both the -i and -o flag can take in gzipped files (see examples in "Sample Commands"). I'm also happy to announce that BetaScan format can now be generated from several file formats, including plink and vcf, by the toolkit [glactools](https://grenaud.github.io/glactools/). Thank you to Gabriel Renaud for these updates!

7/5/17: Added some more checks for valid parameter choice and input file format. Also, slightly modified behavior of script when determining SNPs in the current window, so that core SNP is excluded from window size when the window size is even. This means that the window will now be symmetric around the core SNP, whether an odd or even window size parameter is given. Also, made a tweak so that performance should be slightly quicker than the older version.

5/4/17: Beta now can take in variable sample sizes for each SNP. In other words, not all frequencies have to be calculated using the same number of individuals. Because of this, the input file format has been updated.

## Getting Started
BetaScan is a command line program implemented in python.

### Prerequisites 
* [Python 2](https://www.python.org/downloads/) -Language code is written in
* [Numpy](http://www.numpy.org/) -Python package used for arrays

## Basic Usage

### Input File Format
BetaScan takes in a tab separated file with three columns. The first column contains the coordinate of each variant, and the second contains the frequency of the **derived** allele (note: this is opposite of the BALLET software), in number of haploid individuals, of the variant. However, in practice, for folded Beta only, it doesn't matter if the derived, ancestral, or already folded allele frequency is used in the second column, as BetaScan will fold the frequency anyway. The third column contains the sample size, in number of haploid individuals, that were used to calculate the frequency of that variant. The file should be sorted by position (the unix command sort -g will do this for you). If you are using the Beta2 statistic, substitutions should be coded as SNPs of frequency equal to the sample size.  The scan should be run on each chromosome separately. An example of a sample file is below:

```
14  2 99  
15 99 99
25  1 100  
47  99  100
48  82  95
98 100 100
103 10  100
245 93  96
```
### Parameters 
* -i: Path of input file
* -w: Total window size (default: 1000, corresponding to 500 bp on either side of the center (i.e. core) SNP)
* -p: Value of p (default: 2)
* -m: Minimum folded frequency of core SNP, exclusive, can range from 0 to .5 (default: 0)
* -fold: Use folded version (default: false)
* -o: Output file name (default: print to screen)
* -B2: Calculate the  Beta2 statistic.
* -DivTime: An estimate of the divergence time between your species (only used with the -B2 flag). 
* -SigTest: Instead of just return the Beta value, also return the standardized Beta value. This is equal to the Beta value divided by the square root of the variance. Note that this calculation is slower than calculating just the Beta values, however we do have plans to speed it up.
* -theta: An estimate of the genome-wide mutation rate (only used with the -SigTest option). 

### Explanation of parameters
- w: The sliding window size to use. The expected size of the window in which most of the signal of balancing selection will be found is exponentially distributed with rate rho\*T, where rho is the individual recombination rate, and T is an estimate of the age of the selection in generations. We found that somewhere around the 95th percentile on either side of the core site works well. To calculate the total window size in R, you can use the command 2*qexp(.95, rho*T). Of course, in practice you don't know T, but a ballpark estimate should do.

* -m: In theory, Beta has good performance down to very low and high frequencies. However, the chance of seeing an allele at very low or high equilibrium frequency that is under long-term balancing selection is very low. This is because genetic drift is expected to drift the allele out of the population (see Ewens & Thomson 1970 or Carr & Nassar 1970). We see this phenomenon for variants at folded frequency ~15% or less when we simulate overdominance with human parameters. In addition, poor variant calling can cause false variant calls at low allele frequencies. 

  For this reason, we give the user the option to exclude all variants from consideration at or below a certain folded frequency. These variants will still be used when calculating the Beta score for other SNPs, but their score will not be calculated by the program. On top of decreasing false positives, this has the added benefit of speeding up run-time. Although by default the minimum frequency is zero, we recommend that the user be cautious with interpreting any signal of balancing selection, whether detected using Beta or another method, if the core SNP is at a very extreme frequency.

* -fold: The default version of Beta takes into account the frequency of each variant. However, if ancestral state cannot be confidently called, perhaps due to there being no suitable outgroup, the folded version of Beta should be used. The formulation for this statistic can be found in the supplement of our 2017 paper.

* -B2: This calculates Beta2, which is the new version detailed in our [bioRxiv preprint](https://www.biorxiv.org/content/early/2018/12/17/497255). It has higher power than the Beta1 statistics, but requires substitution data and an estimate of divergence time.

* -DivTime:  An estimate of the divergence time between your two species. This can either be obtained from prior demographic inference of your species, or using the BALLET software. In practice, this value affects the relative ranking of loci very little, but will affect the standardized statistic.  To convert from generations (g) to coalescent units (c), the formula is g=c*Ne*2 where Ne is the effective population size.

* -theta: The estimated mutation rate. It's equal to 2*l*N_e*u, where u is the locus neutral mutation rate, Ne is the effective population size and l is the ploidy. This can be based on prior estimates from your species, or you can estimate it from your own data. We have future plans to automate this calculation, but in the meantime you can calculate Watterson's theta on your data.



### Sample Commands
To run BetaScan on our file SNPFreqs.txt with default parameters:
```
python BetaScan.py -i SNPFreqs.txt
```
To run with a 2000 base pair window, a p parameter value of 50 and using the folded version of Beta and with a specified output file path:
```
python BetaScan.py -i SNPFreqs.txt -w 2000 -p 50 -fold -o /path/chr1.beta.out
```
To run with a 5000 base pair window, a p parameter value of 20 and excluding all core SNPs that are of frequency 10% or less, or of frequency 90% or greater:
```
python BetaScan.py -i SNPFreqs.txt -w 5000 -p 20 -m .1
```
To output a gzipped file:
```
python BetaScan.py -i SNPFreqs.txt | gzip > chr1.beta.out.gz
```
To input a gzipped file:
```
python BetaScan.py -i < (zcat chr1.beta.gz )
```
To calculate standardized Beta2 with an estimated divergence time of 12.5 and an estimated mutation rate of 0.001:
```
python BetaScan.py -B2 -DivTime 12.5 -SigTest -Theta 0.001 -i SNPFreqsWithSubs.txt
```

### Output Format
BetaScan outputs a 2 or 3 column tab delimited file, where the first column contains the coordinate of the core variant, and the second column contains its Beta score. If the -SigTest flag is used, the third column is the value of the standardized statistic. Only core SNPs above your minimum folded frequency threshold, specified by -m, will appear in this file.

*Caution: if variant calls within the specified window size of the core variant are not confident, the value of Beta may not be correct, because of artificially elevated or reduced number of mutations. For this reason, we encourage you to use quality filters.*

## FAQ
#### How many samples are needed to run Beta?

We have found that in human a relatively low number of samples, around 5, is sufficient to detect a fairly large proportion of sites. Maximum power is obtained with a sample size of around 20 (so 10 diploid individuals).

#### How do I choose a value of p?

Results are fairly robust to choice.  However we have found that a value of 2 performs well under a wide array of parameters, so we recommend using 2 (the default as of 12/11/18) unless you have a reason not to.

#### Should I use Beta1, Beta1* or Beta2?
If you only have a folded site frequency spectrum (i.e. you don't know what the ancestral versus derived alleles are) you need to used Beta1*.

If you have an unfolded site frequency spectrum, but don't have substitution information with an outgroup species, use Beta1.

If you have substitution information and a unfolded site frequency spectrum, use Beta2. 

#### I have frequency information I calculated using the --freq command in vcftools. How do I convert the vcf output format to the BetaScan output format?

The toolkit [glactools](https://grenaud.github.io/glactools/) is able to convert between vcfs and BetaScan format and is probably the most robust way to do this. 

Alternatively, you can use the following command in unix:
```
awk -F "\t|:" '(NR>1) && ($6!='0') && ($6!='1') && ($3=='2') {OFS="\t"; print$2,$6*$4,$4}' yourfile.frq
```
This command reformats the .frq file and filters out positions that have more than 2 possible alleles, or are at frequency 0 or 100%. Make sure that you use the fold command if you haven't called ancestral/derived alleles. If you have called them, then this awk script assumes that the derived allele is the first allele listed in the .frq file outputted by vcftoos. Also, please double check that this command outputs the right thing from your .frq file! There could always be variations in the .frq format I don't know about.

#### I ran some simulations using the simulation software SLiM, and want to convert them into BetaScan format. Is there an easy way to do this?

Once again, awk can come to the aid:

```
awk '{OFS="\t"}{if ($1=="Genomes:") exit }(($2=="m3") || ($2=="m1")) && ($8!="100") {print $3,$8,"100"}' SLiMFile.out
```

The first thing to note is that SLiM has more than one output file format, and this awk command only works with the SLiM format, not the ms format. Note, in this example, there's two mutation types simulated: m1 and m3, and both are outputted. You should obviously modify this so it works with your simulation details. This command also assumes a sample size of 100. If this is not your sample size, you should replace "100" with your sample size in quotes.
