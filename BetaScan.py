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
	 SNPi, the index in the array for the current SNP under consideration
	 dataList: the numpy array of all SNP locations & frequencies
	"""

	locSNP = dataList[SNPi,0] #the coordinates of the core SNP
	winStart = locSNP-winSize/2
	firstI= prevStarti + np.searchsorted(dataList[prevStarti:,0],winStart,side='left') #array index of start of window, inclusive
	winEnd = locSNP + winSize/2
	endI = prevEndi - 1 + np.searchsorted(dataList[prevEndi:,0],winEnd,side='right') #array index of end of window, exclusive
 	return  (firstI,endI)


def calc_beta_folded(SNPFreqList, coreFreq, numInd,p):
	"""Calculates the value of d, the similarity measure, times i, the frequency from Siewert et al.
		#SNPFreq: freq of SNP under consideration, ranges from 1 to sample size
		#coreFreq: freq of coresite, ranges from 0 to 1
		#p: the p parameter specificying sharpness of peak
		#numInd: the number of haploid individuals used to calculate frequency of core site
	"""

	if SNPFreqList.size==0:
		return 0
	a1 = np.sum(1./np.arange(1,numInd))
	thetaW = len(SNPFreqList[:,0])/a1
	thetaBNum = np.sum(calcD(SNPFreqList[:,0]/SNPFreqList[:,1],coreFreq,p))

	i = np.arange(1,numInd)
	thetaBDenom = np.sum((1./i)*calcD(i/float(numInd),coreFreq,p))

	thetaB = thetaBNum/thetaBDenom
	return thetaB - thetaW



def calc_beta_unfolded(SNPFreqList, coreFreq, numInd,p):
	"""Calculates the unfolded version of Beta from Siewert and Voight
		For use when the ancestral and derived alleles can be confidently called
	
		Parameters:
		SNPFreqList: a list of frequencies, one for each SNP in the window,
			first column ranges from 1 to number of individuals, second columns is # individuals
		coreFreq: the frequency of the core SNP, must range from 0 to 1, exclusive
		numInd: number of individuals used to calculate the core site frequency
		p: value of parameter p
	"""
	if SNPFreqList.size==0:
		return 0
	a1 = np.sum(1./np.arange(1,numInd))
	thetaW = len(SNPFreqList[:,0])/a1
	thetaBNum = sum(calcD(SNPFreqList[:,0]/SNPFreqList[:,1],coreFreq,p)*SNPFreqList[:,0])
	thetaBDenom = np.sum(calcD(np.arange(1,numInd)/float(numInd),coreFreq,p))
	thetaB = thetaBNum/thetaBDenom
	return thetaB - thetaW


def calc_thetabeta_unfolded(SNPFreqList, coreFreq, numInd,p):
	"""Calculates theta_Beta usign the unfolded SFS
		#SNPFreq: freq of SNP under consideration, ranges from 1 to sample size
		#coreFreq: freq of coresite, ranges from 0 to 1
		#p: the p parameter specificying sharpness of peak
		#numInd: the number of haploid individuals used to calculate frequency of core site
	"""

	if SNPFreqList.size==0:
		return 0
	thetaBNum = np.sum(calcD(SNPFreqList[:,0]/SNPFreqList[:,1],coreFreq,p)*SNPFreqList[:,0])

	thetaBDenom = np.sum(calcD(np.arange(1,numInd)/float(numInd),coreFreq,p))

	thetaB = thetaBNum/thetaBDenom
	return thetaB


def calc_thetabeta_folded(SNPFreqList, coreFreq, numInd,p):
	"""Calculates theta_Beta using the folded SFS
		#SNPFreq: freq of SNP under consideration, ranges from 1 to sample size
		#coreFreq: freq of coresite, ranges from 0 to 1
		#p: the p parameter specificying sharpness of peak
		#numInd: the number of haploid individuals used to calculate frequency of core site
	"""

	if SNPFreqList.size==0:
		return 0
	thetaBNum = np.sum(calcD(SNPFreqList[:,0]/SNPFreqList[:,1],coreFreq,p))
	
	thetaBDenom = np.sum((1./np.arange(1,numInd))*calcD(np.arange(1,numInd)/float(numInd),coreFreq,p))

	thetaB = thetaBNum/thetaBDenom
	return thetaB 



def calc_thetaw_unfolded(SNPFreqList, numInd):
	"""Calculates watterson's theta
	
		Parameters:
		SNPFreqList: a list of frequencies, one for each SNP in the window,
			first column ranges from 1 to number of individuals, second columns is # individuals
		numInd: number of individuals used to calculate the core site frequency
	"""
	if SNPFreqList.size==0:
		return 0

	a1 = np.sum(1./np.arange(1,numInd))

	thetaW = len(SNPFreqList[:,0])/a1
	return thetaW


#Calculates theta_D
def calcThetaD(SNPFreqList,c,n):
	"""
		c: Speciation time in coalescent units
		n: Sample Size
		N: Effective pop size
	"""
	if SNPFreqList.size==0:
		return 0

	S = np.where(SNPFreqList[:,0]==SNPFreqList[:,1])[0].shape[0]
	return S/(c+1./n)


def calcBeta2(SNPFreqList,c,n,coreFreq,p):
	SNPs = SNPFreqList[np.where(SNPFreqList[:,0]!=SNPFreqList[:,1])]
	return calc_thetabeta_unfolded(SNPs,coreFreq,n,p)-calcThetaD(SNPFreqList,c,n)


#Calculates the variance of Theta_S
def calcVarThetaD(c,n,theta):
	i = np.arange(2,n+1)
	x = np.sum(1./(i**2.*(i-1)**2.))
	return (1./(c+1./n))**2.*(theta**2.+c*theta+theta/n+theta**2.*x)


def calcT_B2(SNPFreqList,coreFreq,c,n,p,theta):
	'''
	#coreFreq: freq of SNP under consideration, ranges from 1 to sample size
	#n: sample size of core SNP
	#p: the p parameter specifying sharpness of peak
	#theta: genome-wide estimate of the mutation rate
	'''
	notSubsList_noCore = SNPFreqList[np.where(SNPFreqList[:,0]!=SNPFreqList[:,1])]
	thetaB = calc_thetabeta_unfolded(notSubsList_noCore,coreFreq/n,n,p)
	thetasubs = calcThetaD(SNPFreqList,c,n)
	VarD = calcVarThetaD(c,n,theta)
	VarB = calcVTheta(n,theta,coreFreq,p,False)
	return (thetaB-thetasubs)/math.sqrt(VarD+VarB)


def calcD(freq,x,p):
	"""Calculates the value of d, the similarity measure
		#freq: freq of SNP under consideration, ranges from 0 to 1
		#x: freq of coresite, ranges from 0 to 1
		#p: the p parameter specifying sharpness of peak
	"""
	xf = min(x,1.-x)
	f = np.minimum(freq,1.-freq)
	maxdiff = np.maximum(xf,.5-xf)
	corr = ((maxdiff-np.abs(xf-f))/maxdiff)**p
	return corr 


#Using equation 8 from Achaz 2009
def calcT_unfold(SNPFreqList, coreFreq, SNPn, p, theta):
	"""
		#coreFreq: freq of SNP under consideration, ranges from 1 to sample size
		#SNPn: sample size of core SNP
		#p: the p parameter specifying sharpness of peak
		#theta: genome-wide estimate of the mutation rate
	"""

	x = float(coreFreq)/SNPn

	num = np.sum(SNPFreqList[:,0]/SNPFreqList[:,1]*SNPn*omegai(SNPFreqList[:,0]/SNPFreqList[:,1],SNPn, x,p))
	denom = math.sqrt(an(SNPn,x,p)*theta+ Bn(SNPn,x,p)*theta**2.)

	return num/denom


#Calculates variance of a given estimator of theta, eq 7 from Achaz. 
def calcVTheta(n,theta,coreFreq,p,wattersons):
	"""
		#coreFreq: freq of SNP under consideration, ranges from 1 to sample size
		#n: sample size of core SNP
		#p: the p parameter specifying sharpness of peak
		#theta: genome-wide estimate of the mutation rate
		#Wattersons: whether to calculate wattersons theta instead of 
	"""
	wVector = None
	if wattersons==True:
		wVector = 1./np.arange(1,n)
	else:
		wVector = calcD(np.arange(1,n)/float(n),float(coreFreq)/n,p)
	t1 = np.sum(wVector)**(-2.)
	t2 = theta*np.sum(wVector**2.*np.arange(1,n))

	i = np.arange(1,n)
	s1 = np.sum(wVector**2*i**2*sigma(n,np.column_stack([i,i])))


	coords = np.asarray([(j,i) for i in range(1,n) for j in range(i+1,n)])
	iind = np.asarray([i-1 for i in range(1,n) for j in range(i+1,n)])
	jind = np.asarray([j-1 for i in range(1,n) for j in range(i+1,n)])

	s2 = np.sum(coords[:,0]*coords[:,1]*wVector[iind]*wVector[jind]*sigma(n,coords))

	t3 = theta**2.*(s1+2.*s2)
	return t1*(t2+t3)


def calcVTheta_fold(n,theta,coreFreq,p):
	"""
		#coreFreq: freq of SNP under consideration, ranges from 1 to sample size
		#n: sample size of core SNP
		#p: the p parameter specifying sharpness of peak
		#theta: genome-wide estimate of the mutation rate
		#Wattersons: whether to calculate wattersons theta instead of 
	"""

	wVector = calcD(np.arange(1,n/2+1)/float(n),float(coreFreq)/n,p)
	r = np.arange(1,n/2+1)
	t1 = sum(wVector*(1./r+1./(n-r))*1./(1+(r==n-r)))**-2.
	t2 = sum([wVector[i-1]**2.*(phi(n,i)*theta+rho_p_ii(n,i)*theta**2.) for i in range(1,n/2+1)])
	
	coords = np.asarray([(j,i) for i in range(1,n/2+1) for j in range(1,i)])
	t3 = np.sum(wVector[coords[:,0]-1]*wVector[coords[:,1]-1]*rho_p_ij(n,coords[:,0],coords[:,1])*theta**2.)

	return t1*(t2+2.*t3)


def calcCovFolded(n,theta,coreFreq,p):
	"""
		#coreFreq: freq of SNP under consideration, ranges from 1 to sample size
		#n: sample size of core SNP
		#p: the p parameter specifying sharpness of peak
		#theta: genome-wide estimate of the mutation rate
	"""
	r = np.arange(1,n/2+1)
	wVector = calcD(r/float(n),float(coreFreq)/n,p)
	t1 = 1./sum(wVector*(1./r+1./(n-r))*1./(1.+(r==n-r)))
	t2 = 1./sum((1./r+1./(n-r))*1./(1+(r==n-r)))
	coords = np.asarray([(i,j) for i in range(1,n/2+1) for j in range(1,n/2+1)])
	t3 = np.sum(wVector[coords[:,0]-1]*rho_p_ij(100,coords[:,0],coords[:,1])*theta**2.)

	return t1*t2*t3


def calcVarFoldedBeta(n,theta,coreFreq,p):
	"""
		#coreFreq: freq of SNP under consideration, ranges from 1 to sample size
		#n: sample size of core SNP
		#p: the p parameter specifying sharpness of peak
		#theta: genome-wide estimate of the mutation rate
		#Wattersons: whether to calculate wattersons theta instead of 
	"""
	return calcVTheta_fold(n,theta,coreFreq,p)+calcVTheta(n,theta,coreFreq,p,True)-2.*calcCovFolded(n,theta,coreFreq,p)


def omegai(i,SNPn,x,p):
	"""Calculates 9a
		#i:freq of SNP under consideration, ranges between 0 and 1
		#SNPn: number of chromosomes used to calculate frequency of core SNP
		#x: freq of coresite, ranges from 0 to 1
		#p: the p parameter specifying sharpness of peak
	"""
	n1num = calcD(i,x,p)
	n1denom = np.sum(calcD(np.arange(1.,SNPn)/SNPn,x,p))
	n1 = n1num/n1denom
	n2 = (1./(i*SNPn)) /(np.sum(1./np.arange(1.,SNPn)))
	return n1 - n2


#Eq 12a of Achaz
def phi(n,i):
	#n:sample size
	#i: frequency of SNP, in number of individuals
	return n/((1.+(i==n-i))*i*(n-i))

#eq 12b of Achaz
def rho_p_ii(n,i):
	#n:sample size
	#i: frequency of SNP, in number of individuals
	return (sigma(n,np.column_stack([i,i]))+sigma(n,np.column_stack([n-i,n-i]))+2.*sigma(n,np.column_stack([i,n-i])))/(1.+(i==(n-i)))**2.


#eq 12c of Achaz
def rho_p_ij(n,i,j):
	return (sigma(n,np.column_stack([i,j]))+sigma(n,np.column_stack([i,n-j]))+sigma(n,np.column_stack([n-i,j]))+sigma(n,np.column_stack([n-i,n-j])))/((1.+(i==n-i))*(1.+(j==n-j)))


#Returns alpha_n from Achaz 2009, eq 9b
def an(SNPn,x,p):
	'''
		SNPn: Sample size
		x: frequency, ranges from 0 to 1
		p: value of p parameter
	'''
	i=np.arange(1,SNPn)
	return np.sum(i*omegai(i/float(SNPn),SNPn,x,p)**2.)


#Returns Beta_N from Achaz 2009, eq 9c
def Bn(SNPn,x,p):
	'''
		SNPn: Sample size
		x: frequency, ranges from 0 to 1
		p: value of p parameter
	'''
	
	i = np.arange(1,SNPn)
	n1 = np.sum(i**2.*omegai(i/float(SNPn),SNPn,x,p)**2.*sigma(SNPn,np.column_stack([i,i])))


	coords = np.asarray([(j,i) for i in range(1,SNPn) for j in range(1,i)])
	s2 = np.sum(coords[:,0]*coords[:,1]*omegai(coords[:,0]/float(SNPn),SNPn,x,p)*omegai(coords[:,1]/float(SNPn),SNPn,x,p)*sigma(SNPn,coords))

	n2=2.*s2
	return n1+n2


def calcT_fold(SNPFreqList, coreFreq, SNPn, p, theta):
	"""
		#coreFreq: freq of SNP under consideration, ranges from 1 to sample size
		#SNPn: sample size of core SNP
		#p: the p parameter specifying sharpness of peak
		#theta: genome-wide estimate of the mutation rate
	"""

	x = float(coreFreq)/SNPn
	num = calc_beta_folded(SNPFreqList, x, SNPn,p)
	denom = math.sqrt(calcVarFoldedBeta(SNPn,theta,coreFreq,p))
	return num/denom


#Returns sigma from eq 2 or 3 in Fu 1995
def sigma(n,ij):
	'''
		n: sample size
		ij: 2-d array of integers with 2 cols and no rows
	'''
	np.seterr(all='raise')
	res = np.zeros(ij.shape[0])
	#i must be greater than j
	
	ij[:,0],ij[:,1] = ij.max(axis=1),ij.min(axis=1) #flip coordinates if i is less than j
	ci = np.logical_and(ij[:,0]==ij[:,1], ij[:,0]==n/2)
	if np.any(ci)>0:
		res[ci] = 2.*((Fu_an_vec([n])-Fu_an_vec(ij[ci,0]))/(float(n)-ij[ci,0]))-(1./(ij[ci,0]**2.))

	ci = np.logical_and(ij[:,0]==ij[:,1], ij[:,0]<n/2)
	if np.any(ci)>0:
		res[ci] = Fu_Bn(n,ij[ci,0]+1)

	ci = np.logical_and(ij[:,0]==ij[:,1], ij[:,0]>n/2)
	if np.any(ci)>0:
		res[ci] = Fu_Bn(n,ij[ci,0])-1./(ij[ci,0]**2.)

	ci = np.logical_and(ij[:,0]>ij[:,1], ij[:,0]+ij[:,1]==n)
	if np.any(ci)>0:
		res[ci] = (Fu_an_vec([n])-Fu_an_vec(ij[ci,0]))/(n-ij[ci,0]) + (Fu_an_vec([n])-Fu_an_vec(ij[ci,1]))/(n-ij[ci,1]) - (Fu_Bn(n,ij[ci,0])+Fu_Bn(n,ij[ci,1]+1))/2. - 1./(ij[ci,0]*ij[ci,1])

	ci = np.logical_and(ij[:,0]>ij[:,1], ij[:,0]+ij[:,1]<n)
	if np.any(ci)>0:
		res[ci] = (Fu_Bn(n,ij[ci,0]+1)-Fu_Bn(n,ij[ci,0]))/2.

	ci = np.logical_and(ij[:,0]>ij[:,1], ij[:,0]+ij[:,1]>n)
	if np.any(ci)>0:
		res[ci] = (Fu_Bn(n,ij[ci,1])-Fu_Bn(n,ij[ci,1]+1))/2.-(1./(ij[ci,0]*ij[ci,1]))

	return res



#return a_n from Fu 1995, eq 4
def Fu_an_vec(n):
	a = np.insert(np.cumsum(1./np.arange(1,np.amax(n))),0,0)
	return a[np.asarray(n)-1] #minus one for sum being only to n-1


#returns Beta_n(i) from Fu 1995, eq 5
def Fu_Bn(n,i):
	r = 2.0*n/((n-i+1.)*(n-i)) * (Fu_an_vec([n+1])-Fu_an_vec(i)) - (2./(n-i))
	return r


def main():
	
	#Loads the input parameters given by the user
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", help="Name of input file with all SNPs",type=str,required=True)
	parser.add_argument("-o", help="Output file",type=str,default="/dev/stdout")
	parser.add_argument("-w", help="Maximum Window Size (in bp) to calculate Beta in for a single test SNP",type=int,default=1000)
	parser.add_argument("-p", help="Power to raise difference measure by",type=int,default=2)
	parser.add_argument("-fold", help="Use folded SFS version",action="store_true")
	parser.add_argument("-B2",help="Use the Beta2 statistic. To use this, substiution data with an outgroup is needed.",action="store_true")
	parser.add_argument("-m", help="Minimum folded core SNP frequency, exclusive. Must be between 0 and 0.5.",type=float,default=0)
	parser.add_argument("-SigTest",help="Instead of returning Beta value, return normalized Beta Statistic",default=False,action="store_true")
	parser.add_argument("-theta",help="Estimated genome wide theta value per basepair. Used for calculation of variance. It's equal to 2*l*N_e*u, where u is the locus neutral mutation rate, Ne is the effective population size and l is the ploidy",type=float)
	parser.add_argument("-DivTime",help="Divergence time, in coalescent units, between the two species. Only needed if using B^(2). This can be estimated using the BALLET software, or you can use prior estimates for your species of interest. In practice, this value affects power very little, but will affect the standardized statistic.  To convert from generations (g) to coalescent units (c), the formula is g=c*Ne*2 where Ne is the effective population size.",type=float)
	args = parser.parse_args()
	output = open(args.o,'w')



	#Check for valid file format and parameters
	try:
		SNPs = np.loadtxt(open(args.i,'r'),dtype=float)
	except IOError:
		print sys.exit("Error: Input file cannot be found")
	except:
		print sys.exit("Error: Input file in wrong format")
	if args.m<0 or args.m>.5:
		print sys.exit("Error: Parameter m must be between 0 and 0.5.")
	if args.p<=0:
		print sys.exit("Error: Parameter p must be positive.")
	if len(SNPs.shape)<=1:
		print sys.exit("Error: Because the core SNP is excluded from calculations, there must be at least two SNPs in the input file.")
	if args.SigTest and args.theta==None:
		print sys.exit("Error: In order to normalize Beta statistics, a theta value must be provided using the -theta flag.")
	if args.w<2:
		print sys.exit("Error: Window size must be 2 bp or above. However, you probably want to use a window size much larger than 2.")
	if args.SigTest and args.theta<=0:
		print sys.exit("Error: Theta (population-scaled mutation rate) must be a positive value.")
	if args.p>50:
		print sys.exit("Error: P is too large. Reduce value to prevent python numerical errors. See manual for more information.")
	if args.fold and args.B2:
		print sys.exit("Error: You cannot use both B1* (folded Beta) and B2. B1* is for when you have no outgroup, and B2 is for when you can call substiutions with an outgroup. See manual for guidance about which to use.")
	if args.DivTime>1000:
		print sys.exit("Error: Your divergence time seems very high. Divergence time should be in coalescent units, not generations or years.")
	if args.B2 and not np.any(SNPs[:, 1] == SNPs[:, 2]):
		print sys.exit("Error: You chose to calculate Beta2, but your input file contains no substiutions. If you do not have substiution data, please use Beta1 or Beta1*.")

	if not args.SigTest and args.fold:
		output.write("Position\tBeta1*\n")
	elif args.SigTest and args.fold:
		output.write("Position\tBeta1*\tBeta1*_std\n")
	elif args.SigTest and not args.B2:
		output.write("Position\tBeta1\tBeta1_std\n")
	elif not args.B2:
		output.write("Position\tBeta1\n")
	elif args.B2 and not args.SigTest:
		output.write("Position\tBeta2\n")
	else:
		output.write("Position\tBeta2\tBeta2_std\n")

	if not args.B2 and np.any(SNPs[:, 1] == SNPs[:, 2]):
		SNPs = SNPs[np.where(SNPs[:,1]!=SNPs[:,2])]

	prevStarti = 0
	prevEndi = 0
	for SNPi in range(len(SNPs)):
		loc = int(SNPs[SNPi,0])
		freqCount = float(SNPs[SNPi,1])
		sampleN = int(SNPs[SNPi,2])
		freq = freqCount/sampleN

		if freq<1.0-args.m and freq>args.m and sampleN>3:
			core_loc = SNPs[SNPi,0]
			SNPLocs = SNPs[:,0]
			sI,endI = find_win_indx(prevStarti, prevEndi, SNPi, SNPs, args.w)
			prevStarti = sI
			prevEndi = endI
			B = None
			T = None
			if endI>sI:
				SNPSet = np.take(SNPs,range(sI,SNPi)+range(SNPi+1,endI+1),axis=0)[:,1:]
				if args.fold:
					B = calc_beta_folded(SNPSet,freqCount/sampleN,sampleN,args.p)
				elif not args.fold and not args.B2:
					B = calc_beta_unfolded(SNPSet,freqCount/sampleN,sampleN,args.p)
				elif args.B2:
					B = calcBeta2(SNPSet,args.DivTime,sampleN,freqCount/sampleN,args.p)
				if not args.fold and args.SigTest and not args.B2:
					T = calcT_unfold(SNPSet,freqCount,sampleN,args.p,args.theta*args.w)
				elif args.SigTest and args.fold:
					T = calcT_fold(SNPSet,freqCount,sampleN,args.p,args.theta*args.w)
				elif args.SigTest:
					T = calcT_B2(SNPSet,freqCount,args.DivTime,sampleN,args.p,args.theta*args.w)
			if not args.SigTest:
				output.write(str(loc)+"\t"+str(round(B,6))+"\n")
			else:
				output.write(str(loc)+"\t"+str(round(B,6))+"\t"+str(round(T,6))+"\n")
		elif freq>1.0 or freq<0:
			print sys.exit("Error: Input file contains SNP of invalid frequency on line "+str(SNPi)+".")
		elif freq<1.0-args.m and freq>args.m and sampleN<=3:
			print sys.exit("Error: Sample size must be greater than 3 haploid individuals to make inference, or else theta_beta will always equal theta_watterson's. You may wish to increase the m paramter value to exclude this SNP from being a core SNP.")



if __name__ == "__main__":
    main()
