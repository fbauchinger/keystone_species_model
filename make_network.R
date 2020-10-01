# Franziska Bauchinger
# adapted from "Deciphering microbial interactions and detecting keystone species with co-occurrence networks" (Berry and Widder; 2014)
# https://doi.org/10.3389/fmicb.2014.00219

### compute correlation networks based on output from Lotka-Volterra multispecies model (make_community.R)
### correlation coefficients are compared to a null distribution to assess significance

### this script should be in the same directory as the .RData output from make_community.R
### input parameters, output structure and usage at the bottom


# load libraries
library(permute)
library(igraph)

# compute pairwise correlations
PairCor<-function(Ab, Cortype, alpha, sig, its) {
#Analyze only present species
PRES<-as.vector(apply(Ab,2,sum)!=0)
#Real correlation coefficients of species
CorCoef<-cor(Ab[,PRES], method=Cortype)
# evaluate null distribution by shuffling sites
TMP1<-Ab[,PRES]
TMP3<-array(dim=c(ncol(TMP1), ncol(TMP1), its))
	for (h in 1:its) {
	TMP2<-matrix(nrow=nrow(TMP1), ncol=ncol(TMP1))
		for (i in 1:ncol(TMP2)) TMP2[,i]<-TMP1[shuffle(nrow(TMP1)),i]
		TMP3[,,h]<-cor(TMP2, method=Cortype)
	}
#Generate p values from resampling
pvals<-matrix(nrow=nrow(CorCoef), ncol=ncol(CorCoef))
	for (i in 1:nrow(CorCoef)) {
		for (j in 1:ncol(CorCoef)) {
		posP<-sum(CorCoef[i,j]<sort(TMP3[i,j,]))/its # significance of positive corr.
		negP<-sum(CorCoef[i,j]>sort(TMP3[i,j,]))/its # significance of neg corr.	
		pvals[i,j]<-min(posP, negP)
		}
	} 
# Correct p values for multiple comparison
SigPvals<-p.adjust(pvals[lower.tri(pvals)], method="fdr")
CorCoefAdj<-CorCoef[lower.tri(CorCoef)]

TMP5<-lower.tri(CorCoef)
TMP5[lower.tri(TMP5)]<-SigPvals
TMP6<-TMP5
TMP6[lower.tri(TMP6)]<-SigPvals
TMP6<-t(TMP6)
TMP5[upper.tri(TMP5)]<-TMP6[upper.tri(TMP6)]

CorSig<-CorCoef
CorSig[TMP5>=sig]<-0
diag(CorSig)<-0

# Prepare data for network graphs
expanded<-cbind(expand.grid(a=rownames(CorSig), b=rownames(CorSig)),expand.grid(CorSig), expand.grid(alpha[PRES, PRES]), expand.grid(TMP5))
ex<-subset(expanded, (expanded[,3]!=0))
sign<-vector(length=nrow(ex))
sign[ex[,3]>=0]<-"pos"
sign[ex[,3]<0]<-"neg"
AbsCor<-abs(ex[,3])
ex<-cbind(ex, AbsCor,sign)
colnames(ex)<-c("target", "source", "strength", "alpha","Pvals","strength_abs", "sign") #What? What's alpha?
list(network=ex, IntMat=alpha[PRES, PRES], sigCC= CorSig)
}

# calculate relative abundances from absolute abundance profiles
rel.abundance<-function(Ab){
  Ab_temp<-Ab
  for (r in 1:nrow(Ab)){
    if (sum(Ab[r,]) > 0){
      Ab_temp[r,]<-Ab[r,]/sum(Ab[r,])
   } else Ab_temp[r,]<-Ab[r,]
  }
  return(Ab_temp)
}




######################################################################################################################################
## Main program ##
# load arguments from command line

args=(commandArgs(TRUE))
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }

# read in output from make_community.R
InComm<-list.files(pattern="comm.RData")

# loop through output, if multiple runs were performed
for (i in 1:length(InComm)) {
  load(InComm[i])
  basename=paste( substr(InComm[i],1,(nchar(InComm[i])-11)),"_c_",Cortype,"_its_",its,sep="")
  Ab<-rel.abundance(runA$Ab)
  Ab[is.finite(Ab)==FALSE]<-0
  runB<-PairCor(alpha=runA$alpha, Ab=Ab, Cortype= Cortype, its=its,sig=0.05)

  save(runB, file=paste(basename,"_network.RData",sep=""))
  
}


############################################################################
### INPUT PARAMETERS

# Cortype = correlation type; one of "spearman", "pearson" or "kendall", as specified in the documentation of stats::cor() 
# its = number of permutations performed to calculate null distribution 

## loads output of make_community.R (yourfilename_comm.RData); 
## this should therefore be in the same folder

############################################################################
### OUTPUT

# runB

## $network:
### $target = target species (note that correlation coefficients are not directional, despite this naming convention)
### $source = source species (note that correlation coefficients are not directional, despite this naming convention)
### $strength = correlation coefficient between target and source species
### $alpha = interaction strength between target and source species, as used in LVM
### $Pvals = p-value of correlation coefficient, computed by comparison to a null distribution
### $strength_abs = absolute correlation coefficient between target and source species
### $sign = sign of the correlation coefficient between target and source species

## $IntMat = interaction matrix (N x N matrix; N = number of species simulated)

## $sigCC = correlation coefficients between species (N x N matrix, N = number of species simulated; correlation type according to input parameter) 


############################################################################
# USAGE:
# R CMD BATCH '--args Cortype="spear" its=1000' make_network.R make_network_log.txt

