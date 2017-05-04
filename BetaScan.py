import sys
import numpy as np
from StringIO import StringIO
import argparse
import math
import os

def find_win_indx(prevStarti, prevEndi, SNPi, dataList, winSize):
	"""Takes in the previous indices of the starting and end of the window,
	 then returns the appropriate starting and ending index for the next SNP

	 Parameters
	 prevStarti: the starting index in the array of SNP for the previous core SNP's window, inclusive
	 prevEndi: the ending index in the array for the previous SNP's window, inclusive
	 SNPi, the index in the array for the current SNp under consideration
	 dataList: the numpy array of all SNP locations & frequencies
	"""

	locSNP = dataList[SNPi,0] #the coordinates of the core SNP
	winStart = locSNP-winSize/2
	firstI= np.searchsorted(dataList[prevStarti:,0],winStart,side='left') #coordinate of start of window, inclusive
	winEnd = locSNP + winSize/2
	endI = np.searchsorted(dataList[prevEndi:,0],winEnd,side='left') #coordinate of end of window, exclusive
 	return  (firstI,endI)


def calc_beta_folded(SNPFreqList, coreFreq, numInd,p):
	"""Calculates the value of d, the similarity measure, times i, the frequency from Siewert et al.
		#SNPFreq: freq of SNP under consideration, ranges from 1 to sample size
		#x: freq of coresite, ranges from 0 to 1
		#p: the p parameter specificying sharpness of peak
	"""
	fun = np.vectorize(calcD_fold,otypes=[np.float])
	a1 = 0.
	for i in range(1,numInd):
		a1 += 1./i

	thetaW = len(SNPFreqList[:,0])/a1
	thetaBNum = sum(fun(SNPFreqList[:,0],SNPFreqList[:,1],coreFreq,p))

	thetaBDenom = 0
	for i in range(1,numInd):
		thetaBDenom += (1./i)*calcD_fold(i,numInd,coreFreq,p)

	thetaB = thetaBNum/thetaBDenom

	return thetaB - thetaW



def calcD_fold(SNPFreq,SNPn,x,p):
	"""Calculates the value of d, the similarity measure, times i, the frequency from Siewert et al.
		#SNPFreq: freq of SNP under consideration, ranges from 1 to sample size
		#x: freq of coresite, ranges from 0 to 1
		#p: the p parameter specificying sharpness of peak
	"""
	freq = float(SNPFreq)/SNPn

	x = min(x,1.-x)
	f = min(freq,1.-freq)
	maxdiff = max(x,.5-x)
	corr = (((maxdiff-abs(x-f))/maxdiff)**p)

	return corr 



def calc_beta_unfolded(SNPFreqList, coreFreq, numInd,p):
	"""Calculates the unfolded version of Beta from Siewert et al.
		For use when the ancestral and derived alleles can be confidently called
	
		Parameters:
		SNPFreqList: a list of frequencies, one for each SNP in the window,
			first column ranges from 1 to number of individuals, second columns is # individuals
		coreFreq: the frequency of the core SNP, must range from 0 to 1, exclusive
	"""
	fun = np.vectorize(calcDf_unfold,otypes=[np.float])

	a1 = 0.
	for i in range(1,numInd):
		a1 += 1./i

	thetaW = len(SNPFreqList[:,0])/a1

	thetaBNum = sum(fun(SNPFreqList[:,0],SNPFreqList[:,1],coreFreq,p))
	thetaBDenom = 0
	for i in range(1,numInd):
		thetaBDenom += (calcDf_unfold(i,numInd,coreFreq,p))/float(i)

	thetaB = thetaBNum/thetaBDenom
	return thetaB - thetaW



def calcDf_unfold(SNPFreq,SNPn,x,p):
	"""Calculates the value of d, the similarity measure, times i, the frequency from Siewert et al.
		#SNPFreq: freq of SNP under consideration, ranges from 1 to sample size
		#x: freq of coresite, ranges from 0 to 1
		#p: the p parameter specificying sharpness of peak
	"""
	freq = float(SNPFreq)/SNPn
	x = min(x,1.-x)
	f = min(freq,1.-freq)
	maxdiff = max(x,.5-x)
	corr = (SNPFreq)*(((maxdiff-abs(x-f))/maxdiff)**p)

	return corr 




def main():
	
	#Loads the input parameters given by the user
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", help="Name of input file with all SNPs",type=str,required=True)
	parser.add_argument("-w", help="Maximum Window Size (in bp) to calculate Beta in for a single test SNP",type=int,default=1000)
	parser.add_argument("-p", help="Power to raise difference measure by",type=int,default=20)
	parser.add_argument("-fold", help="Use folded SFS version",action="store_true")
	parser.add_argument("-m", help="Minimum folded core SNP frequency, exclusive",type=float,default=0)

	args = parser.parse_args()


	output = open("Betas_"+args.i.split("/")[-1],'w')
	try:
		SNPs = np.loadtxt(open(args.i,'r'),dtype=float)
	except IOError:
		print sys.exit("Error: Input file cannot be found")
	except:
		print sys.exit("Error: Input file in wrong format")
	prevStarti = 0
	prevEndi = 0
	for SNPi in range(len(SNPs)):
		loc = int(SNPs[SNPi,0])
		freqCount = float(SNPs[SNPi,1])
		sampleN = int(SNPs[SNPi,2])
		freq = freqCount/sampleN
		if freq<1.0-args.m and freq>args.m:
			core_loc = SNPs[SNPi,0]
			SNPLocs = SNPs[:,0]

			sI,endI = find_win_indx(prevStarti, prevEndi, SNPi, SNPs, args.w)
			B = 0
			if endI>sI:
				SNPSet = np.take(SNPs,range(sI,SNPi)+range(SNPi+1,endI),axis=0)[:,1:]
				if args.fold:
					B = calc_beta_folded(SNPSet,freqCount/sampleN,sampleN,args.p)
				else:
					B = calc_beta_unfolded(SNPSet,freqCount/sampleN,sampleN,args.p)

			output.write(str(loc)+"\t"+str(B)+"\n")
		elif freq>1.0 or freq<0:
			print sys.exit("Error: Input file contains SNP of invalid frequency on line "+str(SNPi)+".")




if __name__ == "__main__":
    main()
