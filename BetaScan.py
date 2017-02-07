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

	firstI = prevStarti
	while dataList[firstI,0]<locSNP-winSize/2:
		firstI += 1

	endI = prevEndi
	while endI + 1<len(dataList[:,0]) and dataList[endI+1,0]<locSNP + winSize/2:
		endI +=1


 	return  (firstI,endI)


def calc_beta_folded(SNPFreqList, coreFreq, numInd,p):
	"""Calculates the folded version of Beta from Siewert et al.
	
		Parameters:
		SNPFreqList: a list of frequencies, one for each SNP in the window
		coreFreq: the frequency of the core SNP, must range from 0 to 1, exclusive
	"""
	fun = np.vectorize(calcD_fold,otypes=[np.float])
	a1 = 0.
	for i in range(1,numInd):
		a1 += 1./i

	thetaW = len(SNPFreqList)/a1
	thetaBNum = sum(fun(SNPFreqList,coreFreq,p))

	thetaBDenom = 0
	for i in range(1,numInd):
		thetaBDenom += (1./i)*calcD_fold(float(i)/numInd,coreFreq,p)

	thetaB = thetaBNum/thetaBDenom
	return thetaB - thetaW



def calcD_fold(SNPFreq,x,p):
	"""Calculates the value of d, the similarity measure, from Siewert et al.
		#SNPFreq: freq of SNP under consideration, ranges from 0 to 1
		#x: freq of coresite, ranges from 0 to 1
		#p: the p parameter specificying sharpness of peak
	"""
	x = min(x,1.-x)
	f = min(SNPFreq,1.-SNPFreq)
	maxdiff = max(x,.5-x)
	corr = (((maxdiff-abs(x-f))/maxdiff)**p)
	return corr 



def calc_beta_unfolded(SNPFreqList, coreFreq, numInd,p):
	"""Calculates the unfolded version of Beta from Siewert et al.
		For use when the ancestral and derived alleles can be confidently called
	
		Parameters:
		SNPFreqList: a list of frequencies, one for each SNP in the window
		coreFreq: the frequency of the core SNP, must range from 0 to 1, exclusive
	"""
	fun = np.vectorize(calcDf_unfold,otypes=[np.float])

	a1 = 0.
	for i in range(1,numInd):
		a1 += 1./i

	thetaW = len(SNPFreqList)/a1
	thetaBNum = sum(fun(SNPFreqList,coreFreq,p,numInd))
	thetaBDenom = 0
	for i in range(1,numInd):
		thetaBDenom += (calcDf_unfold(float(i)/numInd,coreFreq,p,numInd))/float(i)

	thetaB = thetaBNum/thetaBDenom
	
	return thetaB - thetaW



def calcDf_unfold(SNPFreq,x,p,numInd):
	"""Calculates the value of d, the similarity measure, times i, the frequency from Siewert et al.
		#SNPFreq: freq of SNP under consideration, ranges from 0 to 1
		#x: freq of coresite, ranges from 0 to 1
		#p: the p parameter specificying sharpness of peak
	"""
	x = min(x,1.-x)
	f = min(SNPFreq,1.-SNPFreq)
	maxdiff = max(x,.5-x)
	corr = (SNPFreq*numInd)*(((maxdiff-abs(x-f))/maxdiff)**p)
	return corr 





def main():
	
	#Loads the input parameters given by the user
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", help="Name of input file with all SNPs",type=str)
	parser.add_argument("-w", help="Maximum Window Size (in bp) to calculate Beta in for a single test SNP",type=int)
	parser.add_argument("-n", help="Number of Individuals in data set",type=int)
	parser.add_argument("-p", help="Power to raise different measure by",type=int)
	parser.add_argument("-fold", help="Use folded SFS version",action="store_true")
	args = parser.parse_args()


	output = open("Betas_"+args.i.split("/")[-1],'w')

	SNPs = np.loadtxt(open(args.i,'r'),dtype=float)
	prevStarti = 0
	prevEndi = 0
	for SNPi in range(len(SNPs)):
		loc = int(SNPs[SNPi,0])
		freq = float(SNPs[SNPi,1])
		if freq<1.0 and freq>0:
			core_loc = SNPs[SNPi,0]
			SNPLocs = SNPs[:,0]

			sI,endI = find_win_indx(prevStarti, prevEndi, SNPi, SNPs, args.w)
			B = 0
			if endI>sI:
				SNPSet = np.take(SNPs,range(sI,SNPi)+range(SNPi+1,endI+1),axis=0)
				if args.fold:
					B = calc_beta_folded(SNPSet[:,1],SNPs[SNPi,1],args.n,args.p)
				else:
					B = calc_beta_unfolded(SNPSet[:,1],SNPs[SNPi,1],args.n,args.p)

			output.write(str(loc)+"\t"+str(B)+"\n")



if __name__ == "__main__":
    main()
