# BetaScan
Welcome! BetaScan implements the β statistic to detect ancient balancing selection, as described in [Siewert & Voight, 2017](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx209/3988103/Detecting-Long-term-Balancing-Selection-using) and [Siewert & Voight, 2020](https://academic.oup.com/gbe/article/12/2/3873/5721358). For in-depth instructions, please read the [BetaScan wiki](https://github.com/ksiewert/BetaScan/wiki).

Update: 11/8/22: The newest version of BetaScan now uses python3! If you want to continue using the old python2 code, you can use BetaScan_python2.py. Otherwise, use BetaScan.py

## Basic Usage
To run BetaScan on an input file named SNPFreqs.txt with default parameters:
```
python BetaScan.py -i SNPFreqs.txt
```
If you have a folded site frequency spectrum, you must include the -fold flag to calculate β<sup>(1)\*</sup>. If your data includes substitutions with an outgroup, you can use the -B2 flag, which calculates β<sup>(2)</sup>. However, if you use -B2 you must include an estimate of speciation time. See the [Usage page of the wiki](https://github.com/ksiewert/BetaScan/wiki/Basic-Usage) for details on how you can estimate.

If you also want to standardize β by its variance, you can do so using the -std flag. This flag must be accompanied by an estimate of the mutation rate using the -theta flag. Once again, see the [Usage page of the wiki](https://github.com/ksiewert/BetaScan/wiki/Basic-Usage) for details on how you can estimate the mutation rate.

## Questions? Comments?
Any feedback or questions are very welcome. You can e-mail Katie at ksiewert@hsph.harvard.edu or post a github issue. We know that programs written by other people can be difficult to use, so we’ve tried our best to make this program simple and intuitive. That being said, bioinformatics is bioinformatics, and issues will arise, so don’t hesitate to contact us!

## References
The original Beta statistics are described in [Detecting Long-Term Balancing Selection Using Allele Frequency Correlation, MBE 2017](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx209/3988103/Detecting-Long-term-Balancing-Selection-using).

Recent updates to BetaScan, including the β<sup>(2)</sup> statistic and standardization are now published in [BetaScan2: Standardized statistics to detect balancing selection utilizing substitution data, GBE 2020](https://academic.oup.com/gbe/advance-article/doi/10.1093/gbe/evaa013/5721358).

## 1000 Genomes Beta Scores
If you would like the β<sup>(1)</sup> Scores for each population in the 1000 Genomes dataset, they are available [here](http://coruscant.itmat.upenn.edu/data/SiewertEA_Full_BetaScores.tar.gz) (warning: this is a 1.8 GB gzipped file). If you just want to look at the top 1% highest scoring haplotypes in each population, that data is available [here](http://coruscant.itmat.upenn.edu/data/SiewertEA_BetaScores.tar.gz). These scores are based on hg19.

β<sup>(2)</sup> Scores are available for the following populations: [YRI](http://coruscant.itmat.upenn.edu/data/YRIStdB2.tar.gz), [CEU](http://coruscant.itmat.upenn.edu/data/CEUStdB2.tar.gz) and [CHB](http://coruscant.itmat.upenn.edu/data/CHBStdB2.tar.gz). And here's the [ReadMe.txt](http://coruscant.itmat.upenn.edu/data/README_B2stdscores.txt). These scores are also based on hg19.

