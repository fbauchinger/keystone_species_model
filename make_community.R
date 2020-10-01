# Franziska Bauchinger
# adapted from "Deciphering microbial interactions and detecting keystone species with co-occurrence networks" (Berry and Widder; 2014)
# https://doi.org/10.3389/fmicb.2014.00219

### Lotka-Volterra multispecies model
### input parameters, output structure and usage at the bottom

# load libraries
library(deSolve)
library(permute)
library(stats)
library(igraph)
library(seqtime)



## Model functions
# Multispecies Lotka-Volterra model
lvm <- function(t,x,parms){
	with(as.list(parms,c(x)),{
		x[x<0.001]<-0 # If abundance is less than 0.001 round to zero
	 dx <-((r*x)*(1-(a %*% x)/k))
	     list(dx)
	})
}

# Numerical integration of mLVM
n.integrate <- function(time=time,init.x=init.x,model=model, parms=parms){
        t.out <- seq(time$start,time$end,length=time$steps)
	as.data.frame(lsoda(init.x,t.out,lvm,parms))
}

# Generate random values for coefficients
generate.parameters <- function(C,IntType,n,mode,K,sites,neg,adjIA){
	length.index <- length(n)
	
	if (IntType == "klemm1"){
	  IntType<-"klemm"
	  clique.size<-2
	}
	
	if (IntType == "klemm2"){
	  IntType<-"klemm"
	  clique.size<-8
	}
	
  if (mode== 1){
    
    params<-mode1(C=C,IntType=IntType,n=n,neg=neg,clique.size=clique.size,K=K)
    
  }
	
	if (mode== 2){
	  
	  params<-mode2(C=C,IntType=IntType,n=n,neg=neg,clique.size=clique.size,K=K)
	  
	}
	
	if (mode== 3){
	  
	  params<-mode3(C=C,IntType=IntType,n=n,neg=neg,clique.size=clique.size,K=K)
	  
	}
	
	alpha<-params$alpha
	r<-params$r
	cc<-params$cc
	
  e<-r
  
	if (adjIA == TRUE){
	  
	  outIA<-apply(alpha != 0,2,sum)
	  inIA<-apply(alpha != 0,1,sum)
	  numIA<-outIA+inIA
	  top<-sort(numIA,decreasing=TRUE)[1:6]
	  pos<-match(top,numIA)
	  pos<-sample(pos,6,replace=FALSE)
	  in.n<-pos[1:3]
	  out.n<-pos[4:6]
	  
	  for (i in out.n){
	    # make out-going interactions stronger
	    z<-length(which(alpha[,i] < 0))
	    if (z!=0) {alpha[alpha[,i]<0,i]<-runif(n=z,min=-0.9,max=-0.7)}
	    z<-length(which(alpha[,i] > 0))
	    if (z!=0) {alpha[alpha[,i]>0,i]<-runif(n=z,min=0.7,max=0.9)}
	    # make in-going interactions weaker
	    z<-length(which(alpha[i,] < 0))
	    if (z!=0) {alpha[i,alpha[i,]<0]<-runif(n=z,min=-0.3,max=-0.1)}
	    z<-length(which(alpha[i,] > 0))
	    if (z!=0) {alpha[i,alpha[i,]>0]<-runif(n=z,min=0.1,max=0.3)}
	  }
	  
	  for (i in in.n){
	    # make out-going interactions weaker
	    z<-length(which(alpha[,i] < 0))
	    if (z!=0) {alpha[alpha[,i]<0,i]<-runif(n=z,min=-0.3,max=-0.1)}
	    z<-length(which(alpha[,i] > 0))
	    if (z!=0) {alpha[alpha[,i]>0,i]<-runif(n=z,min=0.1,max=0.3)}
	    # make in-going interactions stronger
	    z<-length(which(alpha[i,] < 0))
	    if (z!=0) {alpha[i,alpha[i,]<0]<-runif(n=z,min=-0.9,max=-0.7)}
	    z<-length(which(alpha[i,] > 0))
	    if (z!=0) {alpha[i,alpha[i,]>0]<-runif(n=z,min=0.7,max=0.9)}
	  }
	  
	}

	# make sure that diagonal elements are 1
	diag(alpha) <- 1

  list(aMatrix=alpha,growth=r,capacity=cc,extinction=e)
}


#### mode 1
mode1<-function(C,IntType,n,neg,clique.size){
  
  length.index <- length(n)
  
  # create interaction matrix, for details see seqtime documentation
  alpha <- as.matrix(generateA(N=length.index,type=IntType,c=C,pep=(100-neg),max.strength=0.99,min.strength=-0.99,clique.size=clique.size))

  r <- numeric(length.index)
  cc<- numeric(length.index)
  
  # generate growth rates/immigration probabilities
  r[n]=as.numeric(rep(0.5,length.index))
  r[r>=1]<-0.99
  r[r<=0]<-0.01
  r<-r[shuffle(r)]
  
  # generate carrying capacities for species with indices in index vector using beta distribution
  cc[n]<-as.numeric(rep(K,length.index))
  cc[cc<0]<-0 #remove negative values
  cc<-cc[shuffle(cc)]
  
  return(list(alpha=alpha,r=r,cc=cc))
  
}


#### mode 2
mode2<-function(C,IntType,n,neg,clique.size){
  
  length.index <- length(n)
  
  # create interaction matrix, for details see seqtime documentation
  alpha <- as.matrix(generateA(N=length.index,type=IntType,c=C,pep=(100-neg),max.strength=0.99,min.strength=-0.99,clique.size=clique.size))
  
  r <- numeric(length.index)
  cc<- numeric(length.index)
  
  # generate growth rates/immigration probabilities
  r[n]=as.numeric(rexp(length.index,0.5))
  r[r>=1]<-0.99
  r[r<=0]<-0.01
  r<-r[shuffle(r)]
  
  # generate carrying capacities for species with indices in index vector using beta distribution
  cc[n]<-as.numeric(rep(K,length.index))
  cc[cc<0]<-0 #remove negative values
  cc<-cc[shuffle(cc)]
  
  return(list(alpha=alpha,r=r,cc=cc))
  
}

#### mode 3
mode3<-function(C,IntType,n,neg,clique.size){
  
  length.index <- length(n)
  
  # create interaction matrix, for details see seqtime documentation
  alpha <- as.matrix(generateA(N=length.index,type=IntType,c=C,pep=(100-neg),max.strength=0.99,min.strength=-0.99,clique.size=clique.size))
  
  r <- numeric(length.index)
  cc<- numeric(length.index)
  
  # generate growth rates/immigration probabilities
  r[n]=as.numeric(rexp(length.index,0.5))
  r[r>=1]<-0.99
  r[r<=0]<-0.01
  r<-r[shuffle(r)]
  
  # generate carrying capacities for species with indices in index vector using beta distribution
  #count<-sample(3000:5000,1)
  cc[n]<-as.numeric(round(generateAbundances(N=length.index,count=K,mode=5)))
  cc[cc<0]<-0 #remove negative values
  cc<-cc[shuffle(cc)]
  
  return(list(alpha=alpha,r=r,cc=cc))
  
}




# Generate species abundance table for multiple sampling sites with glv
LVMits<-function(time,parameters,sharedSp, sites,nMeta,K,pInt) {
	r <- parameters$growth
	a <- parameters$aMatrix
	cc <- parameters$capacity
	#Set carrying capacity to essentially zero to avoid erros with integrator
	cc[cc==0]<-0.00001
	# Subsample metacommunity and run multiple iterations
	outt<-matrix(nrow= sites, ncol=nMeta)
	rownames(outt)<-paste("Site",seq(1:sites), sep=" ")
	colnames(outt)<-paste("Species",seq(1:nMeta), sep=" ")
	init<-matrix(nrow=sites, ncol=nMeta)
		for (i in 1:sites) {
		if (i<(sites/2)) {ssInd<-c(rbinom(round(nMeta/2), 1, sharedSp*(pInt+1)),rbinom((nMeta-round(nMeta/2)), 1, sharedSp*(1-pInt)))
			} else partition<-ssInd<-c(rbinom(round(nMeta/2), 1, sharedSp*(1-pInt)),rbinom((nMeta-round(nMeta/2)), 1, sharedSp*(pInt+1)))
		parms<-list(r=r, a=a, k=cc)
		init.x<-generateAbundances(nMeta,mode=2,count=K)
		init.x<-init.x[shuffle(init.x)]*ssInd
		init[i,]<-init.x
		# run actual lvm model
		out <-n.integrate(time=time,init.x=init.x,model=lvm, parms=parms)
		outt[i,]<-as.numeric(out[nrow(out),2:ncol(out)])
		}
list(Ab=outt, alpha=a, growth=r, carcap=cc, init=init, time=time)
}




######################################################################################################################################
## Main program ##

#load arguments from command line
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

nMeta<-round(siteRich/sharedSp,digits=0)

for (i in 1:length(IntType)) { # loop through multiple network times
  for (j in 1:length(C)){      # loop through multiple connectivities
    
    # generate parameters for lvm model, according to parameters provided
    parameters<-generate.parameters(C=C[j], IntType=IntType[i], n=c(1:nMeta), mode=mode, K=K, sites=sites, neg=neg, adjIA=adjIA)
    # run lvm model
    runA<-LVMits(time=time,parameters=parameters, sharedSp=sharedSp, sites=sites, nMeta=nMeta, K=K, pInt=pInt)
    # save output in .RData format, output name lists parameters specified
    save(runA, file=paste("Model_glv_C_",C[j],"_neg_",neg,"_mode_",mode,"_adjIA_",adjIA,"_R_",siteRich,"_ss_",sharedSp,"_s_",sites,"_IntType_",IntType[i],"_comm.RData",sep=""))
  }
}


############################################################################
### INPUT PARAMETERS

# time = list with "start", "end" and "steps"
# C = connectivity of the interaction network (overall connectivity, independent of structure), 
#     can also be a vector of values (program will loop through)
#     in a highly structured network (klemm1 and klemm2) lower connectivities can already lead to abundance explosions
# siteRich = average number of species at a given site
# sharedSp = proportion of species at a site that are shared with at least one other site
# sites = number of sites simulated (each site is a separate model run)
# IntType = type of network: "random" for uniformly distributed interactions
#                            "klemm1" for a klemm-eguiluz network with clique.size = 2 (see seqtime documentation)
#                            "klemm2" for a klemm-eguiluz network with clique.size = 8 (see seqtime documentation)
#           can also be a vector of network types (programm will loop through)
# mode = gives the structure of the model parameters
#        "mode1": growth rate at 0.5 for all species, carrying capacity at 50 for all species, network according to IntType with IA-strength [-0.99;0.99]
#        "mode2": carrying capacities and IA-Network as mode1, growth rates drawn from an exponential distribution with values [0;1]
#        "mode3": IA-Network as mode1, growth rates as mode2, carrying capacities drawn from a broken stick distribution (function in vegan) with count between 3000 and 5000
# K = depending on the mode, either carrying capacity for all species or maximum carrying capacity
# pInt = partition analysis. Split sites in two equal sized groups and partition species between them with an intensity of 0 (none) to 1 (perfect exclusion). (PInt+1)*SharedSp must be less than or equal to 1!!!!
# neg = negative edge percentage of the interaction network (NOTE: negative edges are positive interactions!!, more negative interactions can prevent abundance explosions to some extent)
# adjIA = if TRUE, the nodes with the most edges are adjusted so they have either particularly strong out- or in-going interactions (keeping the overall interactions the same)


############################################################################
### OUTPUT
# runA
## $Ab = abundance table; sites in rows, species in columns
## $alpha = interaction matrix (N x N; N = number of species simulated)
## $growth = growth rates per species (as specified by "mode" in input parameters)
## $carcap = carrying capacity (as specified by "mode" in input parameters)
## $init = initial abundances for each species and each simulated site (sampled from a uniform distribution, for details see documentation of seqtime::generateAbundances)
## $time = list with "start"-timepoint, "end"-timepoint and "steps" (as specified with input parameters)



############################################################################
# USAGE:
# R CMD BATCH '--args time=list(start=0,end=200,steps=1000) C=c(0.02,0.03) siteRich=10 sharedSp=0.8 sites=10 IntType=c("klemm1","random") mode=1 K=10 pInt=0 neg=40 adjIA=FALSE' make_community.R make_community_log.txt


