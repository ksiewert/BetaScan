# BetaScan
BetaScan implements the Beta statistic to detect ancient balancing selection. Our paper describing this statistic and its application is available [here](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx209/3988103/Detecting-Long-term-Balancing-Selection-using). BetaScan takes in a file of variant positions and frequencies, and calculates Beta in a sliding window approach. It allows the user to choose appropriate parameter choices, and outputs the value of Beta for each variant.

Any feedback or questions are very welcome. You can e-mail Katie at ksiewert@pennmedicine.upenn.edu or post a github issue. We know that programs written by other people can be difficult to use, so we’ve tried our best to make this program simple and intuitive. That being said, bioinformatics is bioinformatics, and issues will arise, so don’t hesitate to contact us!

# 1000 Genomes Beta Scores
If you would like the Beta Scores for each population in the 1000 Genomes dataset, they are available [here](http://coruscant.itmat.upenn.edu/data/SiewertEA_Full_BetaScores.tar.gz) (warning: this is a 1.8 GB gzipped file). If you just want to look at the top scoring haplotypes in each population, that data is available [here](http://coruscant.itmat.upenn.edu/data/SiewertEA_BetaScores.tar.gz).

# Recent Updates
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
BetaScan takes in a tab separated file with three columns. The first column contains the coordinate of each variant, and the second contains the frequency of the **derived** allele (note: this is opposite of the BALLET software), in number of haploid individuals, of the variant. However, in practice, for folded Beta only, it doesn't matter if the derived, ancestral, or already folded allele frequency is used in the second column, as BetaScan will fold the frequency anyway. The third column contains the sample size, in number of haploid individuals, that were used to calculate the frequency of that variant. The file should be sorted by position (the unix command sort -g will do this for you). Variants with frequencies of exactly 0% or 100% should not be included.  The scan should be run on each chromosome separately. An example of a sample file is below:

```
14  2 99  
25  1 100  
47  99  100
48  82  95
103 10  100
245 93  96
```
### Parameters 
* -i: Path of input file
* -w: Total window size (default: 1000, corresponding to 500 bp on either side of the center (i.e. core) SNP)
* -p: Value of p (default: 20)
* -m: Minimum folded frequency of core SNP, exclusive, can range from 0 to .5 (default: 0)
* -fold: Use folded version (default: false)
* -out: Output file name (default: print to screen)

### Explanation of parameters
* -m: In theory, Beta has good performance down to very low and high frequencies. However, the chance of seeing an allele at very low or high equilibrium frequency that is under long-term balancing selection is very low. This is because genetic drift is expected to drift the allele out of the population (see Ewens & Thomson 1970 or Carr & Nassar 1970). We see this phenomenon for variants at folded frequency ~15% or less when we simulate overdominance with human parameters. In addition, poor variant calling can cause false variant calls at low allele frequencies. 

  For this reason, we give the user the option to exclude all variants from consideration at or below a certain folded frequency. These variants will still be used when calculating the Beta score for other SNPs, but their score will not be calculated by the program. On top of decreasing false positives, this has the added benefit of speeding up run-time. Although by default the minimum frequency is zero, we recommend that the user be cautious with interpreting any signal of balancing selection, whether detected using Beta or another method, if the core SNP is at a very extreme frequency.

* -fold: The default version of Beta takes into account the frequency of each variant. However, if ancestral state cannot be confidently called, perhaps due to there being no suitable outgroup, the folded version of Beta should be used. The formulation for this statistic can be found in the supplement of our paper.

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

### Output Format
BetaScan outputs a 2-column tab delimited file, where the first column contains the coordinate of the core variant, and the second column contains its Beta score. Only core variants above your minimum folded frequency threshold, specified by -m, will appear in this file.

*Caution: if variant calls within the specified window size of the core variant are not confident, the value of Beta may not be correct, because of artificially elevated or reduced number of mutations. For this reason, we encourage you to use quality filters.*

## FAQ
#### How many samples are needed to run Beta?

We have found that in human a relatively low number of samples, around 5, is sufficient to detect a fairly large proportion of sites. Maximum power is obtained with a sample size of around 20 (so 10 diploid individuals).

#### How do I choose a value of p?

Although there is no definite best choice of p, and results are fairly robust to choice, we do have some guidelines:
If frequency calls can be called confidently, a value around 20 is sufficient.
If frequency calls are not confident, a smaller value of p (around 2) should be used. 

#### What window size should I use?

See the supplement of our paper on BioRxiv for a rough derivation of maximum window size, based on the estimated recombination rate.

#### Should I use the folded or unfolded version of Beta?

If you have accurate ancestral calls, then we recommend you use the unfolded version, because it can detect balanced haplotypes at more extreme frequencies. If you're not confident in the ancestral calls, then the folded version should be used, because it has identical power throughout most of the site frequency spectrum. In practice, it's not a bad idea to do both. First using unfolded Beta for your general scan, and then using the folded version to double check that any intermediate-frequency top unfolded scan hits are not an artifact of ancestral allele misidentification.

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
