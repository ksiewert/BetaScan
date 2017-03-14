# BetaScan
Beta scan implements the Beta statistic to detect ancient balancing selection. A copy of the paper describing this statistic and its application is available [here on BioRxiv](http://biorxiv.org/content/early/2017/03/01/112870). Beta scan takes in a file of variant positions and frequencies, and calculates Beta in a sliding window approach. It allows the user to choose appropriate parameter choices, and outputs the value of Beta for each variant.

Any feedback or questions are very welcome. You can e-mail Katie Siewert at ksiewert@upenn.edu. We know that programs written by other people can be difficult to use, so we’ve tried our best to make this program simple and intuitive. That being said, bioinformatics is bioinformatics, and issues will arise, so don’t hesitate to contact us!

Caution: This program is in its early stages. I will continue to test it and make changes throughout the upcoming weeks.

## Getting Started
Beta scan is a command line program implemented in python.

### Prerequisites 
* [Python 2](https://www.python.org/downloads/) -Language code is written in
* [Numpy] (http://www.numpy.org/) -Python package used for arrays

## Basic Usage

### Input File Format
Beta Scan takes in a tab separated file with two columns. The first column contains the coordinate of each variant, and the second contains the frequency, between 0 and 1, exclusive, of the variant. The file should be sorted by position (the unix command sort -g will do this for you). Variants with frequencies of exactly 0 or 1 should not be included. The scan should be run on each chromosome separately. An example of a sample file is below:

```
1	.2  
25	.01  
47	.16  
48	.84  
103	.10  
245	.99  
```
### Parameters 
* -i: Path of input file
* -n: Sample size
* -w: Window size (default: 1000)
* -p: Value of p (default: 20)
* -f: Use folded version (default: false)
* -m: Minimum folded frequency of core SNP, exclusive (default: 0)

### Explanation of parameters
-m: In theory Beta has good performance down to very low and high frequencies. However, the chance of seeing an allele at very low or high equilibrium frequency that is under long-term balancing selection is very low. This is because genetic drift is expected to drift the allele out of the population (see Ewens & Thomson 1970 or Carr & Nassar 1970). We see this phenomena for variants at folded frequency ~15% or less when we simulate overdominance with human parameters. In addition, poor variant calling can cause false variant calls at low allele frequencies. 

For this reason, we give the user the option to exclude all variants from consideration at or below a certain folded frequency. These variants will still be used when calculating the Beta score for other SNPs, but their score will not be calculated by the program. On top of decreasing false positives, this has the added benefit of speeding up run-time. Although by default it is zero, we recommend that the user be cautious with interpreting any signal of balancing selection, whether detected using Beta or another method, where the core SNP is at a very extreme frequency.

-fold: The default version of Beta takes into account the frequency of each variant. However, if ancestral state cannot be confidently called, perhaps due to there being no suitable outgroup, the folded version of Beta should be used. The formulation for this statistic can be found in the supplement of our paper.

### Sample Commands
To run Beta Scan on our file SNPFreqs.txt with default parameters:
```
python BScan.py -i SNPFreqs.txt -n 100
```
To run with a 2000 basepair window, a p parameter value of 50 and using the folded version of Beta:
```
python BScan.py -i SNPFreqs.txt -n 100 -w 2000 -p 50 -f
```
To run with a 5000 basepair window, a p parameter value of 20 and excuding all core SNPs that are of frequency 10% or less:
```
python BScan.py -i SNPFreqs.txt -n 100 -w 5000 -p 50 -m .1
```

### Output Format
Beta Scan outputs a 2-column tab delimited file, where the first column contains the coordinate of the core variant, and the second column contains its Beta score.

*Caution: if variant calls within the specified window size of the core variant are not confident, the value of Beta may not be correct, because of artificially elevated or reduced number of mutations. For this reason, we encourage you to use quality filters.*

## FAQ
1. How many samples are needed to run Beta?

We have found that in human, a relatively low number of samples, around 5, are sufficient to detect a fairly large proportion of sites. Maximum power is obtained with a sample size of around 20 (so 10 diploid individuals).

2. How do I choose a value of p?

Although there is no definite best choice of p, and results are fairly robust to choice, we do have some guidelines:
If frequency calls can be called confidently, a value around 20 is sufficient
If frequency calls are not confident, a smaller value of p should be used. 

3. What window size should I use?

See the supplement of the future Siewert & Voight, 2017 for a rough derivation of maximum window size, based on the estimated recombination rate.


