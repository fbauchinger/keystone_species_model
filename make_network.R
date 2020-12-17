# Franziska Bauchinger
# adapted from "Deciphering microbial interactions and detecting keystone species with co-occurrence networks" (Berry and Widder; 2014)
# https://doi.org/10.3389/fmicb.2014.00219


### compute correlation networks based on output from Lotka-Volterra multispecies model (make_community.R)
### correlation coefficients are compared to a null distribution to assess significance

### this script should be in the same directory as the .RData output from make_community.R

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

## $correlation.stats = statistics computed from the correlation network sigCC

### network statistics computed are:
#### rel.abun = mean relative abundance of a species across all sites
#### cstr = sum of absolute correlation strengths of a species
#### mean.cstr = mean of absolute correlation strengths of a species
#### cb = betweenness centrality, calculated with igraph::betweenness
#### cc = closeness centrality, calculated with igraph::closeness
#### cc_rel = relative closeness centrality; normalized to the number of edges in the network
#### md = node degree, calculated with graph::degree
#### clust = clustering coefficient or local transitivity, calculated with igraph::transitivity with type=local
#### eigen.centr = eigenvector centralities, calculated with igraph::eigen_centrality
#### alpha.centr = alpha centrality, calculated with igraph:alpha_centrality
#### transitivity = weighted transitivity, calculated with igraph::transitivity with type = "weighted"

############################################################################
# USAGE:
# R CMD BATCH '--args Cortype="spear" its=1000' make_network.R make_network_log.txt



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


make_active_net<-function(alpha, Ab, nNets, verb){
  adj<-as.matrix(ifelse (alpha!=0 ,1,0))
  g<-graph.adjacency(adj, mode="undirected")
  t<-which(alpha!=0 ,arr.ind=T)
  #reorder by row (just as igraph does obviously)
  t<-t[ order(t[,1]),]
  E(g)$color<-ifelse(alpha[t]>0,"green","red")
  vcol<-matrix(ncol=ncol(Ab),nrow=nNets)
  graphs<-list()
  for(i in 1:nNets){
    graphs[[i]]<-g
    V(graphs[[i]])$color<-ifelse(Ab[i,]>0.01,1,0)
    if(verb==T){plot(graphs[[i]])}
    vcol[i,]<-as.numeric(V(graphs[[i]])$color)
  }
  list(graphs=graphs, Vcolor=vcol)
}


centrality.norm<-function(graph,type=c("degree","closeness","betweenness"),centralization=FALSE)
{
  result<-NA
  g<-graph
  cent<-centralization
  if (!is.igraph(g)) {stop("Not a graph object")}
  if (type[[1]] == "degree") {
    if (!cent) result <- degree(gimode="all")/(vcount(g)-1)
    else result <- (sum(max(degree(g,mode="all"))-degree(g,mode="all")))/((vcount(g)-1)*(vcount(g)-2))}
  else if (type[[1]] == "betweenness") {
    temp <- 2*betweenness(g,directed=F)/((vcount(g)-1)*(vcount(g)-2))
    if (!cent) result <- temp
    else  result <- sum(max(temp)-temp)/(vcount(g)-1)}
  else if (type[[1]] == "closeness") {
    if (!cent) result <- closeness(g,mode="all")
    else result <- (2*vcount(g)-3)*(sum(max(closeness(g,mode="all"))-closeness(g,mode="all")))/((vcount(g)-1)*(vcount(g)-2))}
  else {stop("this type is unavailable or mispelled")}
  return(result)
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
  num.sites<-nrow(Ab)
  runB<-PairCor(alpha=runA$alpha, Ab=Ab, Cortype=Cortype, its=its,sig=0.05)

  Corr<-runB$sigCC
  g<-make_active_net(Corr, Ab, num.sites, verb=F)
  g2<-graph.adjacency(Corr,mode = "undirected", weighted = TRUE)
  
  cbtmp<-centrality.norm(g$graphs[[1]],type="betweenness",centralization=F)
  cctmp<-centrality.norm(g$graphs[[1]],type="closeness",centralization=F)
  cc_reltmp<-(centrality.norm(g$graphs[[1]],type="closeness",centralization=F))/max(centrality.norm(g$graphs[[1]],type="closeness",centralization=F))
  mdtmp<-(degree(g$graphs[[1]])-2)
  clusttmp<-transitivity(g$graphs[[1]],type="local")
  eigencentr<-eigen_centrality(g$graphs[[1]],directed=TRUE,weights=NA)$vector
  alphacentr<-alpha_centrality(g2)
  trans<-transitivity(g2,type = "weighted")
  cstr<-rowSums(abs(Corr))
  mean.cstr<-rowMeans(abs(Corr))
  
  correlation.stats<-data.frame(OTU=rownames(Corr), rel.abun=colMeans(Ab),cstr=cstr,mean.cstr=mean.cstr,cb=cbtmp,cc=cctmp,cc.rel=cc_reltmp,
                                md=mdtmp,clust=clusttmp,eigen.centr=eigencentr,alpha.centr=alphacentr,trans=trans)
  
  runB$correlation.stats<-correlation.stats
  
  save(runB, file=paste(basename,"_network.RData",sep=""))
  
}




