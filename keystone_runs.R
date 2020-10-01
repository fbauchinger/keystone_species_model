# Franziska Bauchinger
# adapted from "Deciphering microbial interactions and detecting keystone species with co-occurrence networks" (Berry and Widder; 2014)
# https://doi.org/10.3389/fmicb.2014.00219

### re-runs multispecies Lotka-Volterra model as specified in make_community.R
### in each run, abundance of one species is set to 0
### the distance between the initial model run, with the species present, and the model run with the species absent, is calculated
### this procedure is repeated multiple times and an average distance for each species is calculated
### this average distance is taken as a proxy for a species keystone potential

### this script should be in the same directory as the .RData output from make_community.R
### input parameters, output structure and usage at the bottom


# load libraries
library(permute)
library(network)
library(igraph)
library(deSolve)
library(vegan)
library(entropy)
library(seqtime)

##iteratively sets every species to zero and returns object similar 
#to run to be pushed into make_active_net with glv
LVMKeys<-function(time, param, init, AB) {
  # Number of species in community
  nMeta <- length(AB)
  # run integration with IC for each species
  outt<-matrix(nrow= (nMeta+1), ncol=nMeta)
  rownames(outt)<-paste("ZeroSpec_", seq(1:(nMeta+1))-1, sep="")
  colnames(outt)<-paste("Species",seq(1:nMeta), sep="")
  outt[1,]<-as.numeric(AB)
  init.x=matrix(nrow= (nMeta+1), ncol=nMeta)
  for( i in 1:(nMeta+1)){init.x[i,]<-as.numeric(init)}
  for( i in 2:(nMeta+1)){init.x[i,i-1]<-0}
  for (i in 2:(nMeta+1)) {
    out <-n.integrate(time=time,init.x=init.x[i,],model=lvm, parms=param)
    outt[i,]<-as.numeric(out[nrow(out),2:ncol(out)])
  }
  list(Ab=outt, init=init.x, growth=param$r, carcap=param$k, alpha=param$a)
}

##########
#make active net takes keyruntmp and calculates the active nodes; returns object
#containing all active graphs and Vertex coloring sheme for each of the graphs $Vcolor
make_active_net<-function(run, nNets, verb){
  adj<-as.matrix(ifelse (run$alpha!=0 ,1,0))
  g<-graph.adjacency(adj, mode="directed")
  t<-which(run$alpha!=0 ,arr.ind=T)
  #reorder by row (just as igraph does obviously)
  t<-t[ order(t[,1]),]
  E(g)$color<-ifelse(run$alpha[t]>0,"green","red")
  vcol<-matrix(ncol=ncol(run$Ab),nrow=nNets)
  graphs<-list()
  for(i in 1:nNets){
    graphs[[i]]<-g
    V(graphs[[i]])$color<-ifelse(run$Ab[i,]>0.01,1,0)
    if(verb==T){plot(graphs[[i]])}
    vcol[i,]<-as.numeric(V(graphs[[i]])$color)
  }
  list(graphs=graphs, Vcolor=vcol)
}
 


# Multispecies Lotka-Volterra model
lvm <- function(t,x,parms){
  with(as.list(parms,c(x)),{
    x[x<0.001]<-0 # If abundance is less than 1 round to zero
    dx <-((r*x)*(1-(a %*% x)/k))
    list(dx)
  })
}


# Numerical integration of MLV model
n.integrate <- function(time=time,init.x=init.x,model=model, parms=parms){
  t.out <- seq(time$start,time$end,length=time$steps)
  as.data.frame(lsoda(init.x,t.out,lvm,parms))
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

keystone.distance<-function(Ab){
  Ab[is.finite(Ab)==FALSE]<-0
  distance_bray<-c()
  distance_canberra<-c()
  for (row in 1:ncol(Ab)) {
    spec<-row+1  
    dist1<-vegdist(Ab[c(1,spec),-row],method="bray")
    distance_bray<-c(distance_bray,dist1)
    dist2<-vegdist(Ab[c(1,spec),-row],method="canberra")
    distance_canberra<-c(distance_canberra,dist2)
  }
  
  ### normalize change in abundance for starting abundance
  for (col in 1:ncol(Ab)){
    Ab[,col]<-Ab[,col]/Ab[1,col]
  }
  Ab[is.finite(Ab)==FALSE]<-0
  ## calculate distances
  distance_norm_bray<-c()
  distance_norm_canberra<-c()
  for (row in 1:ncol(Ab)) {
    spec<-row+1  
    dist1.norm<-vegdist(Ab[c(1,spec),-row],method="bray")
    distance_norm_bray<-c(distance_norm_bray,dist1.norm)
    dist2.norm<-vegdist(Ab[c(1,spec),-row],method="canberra")
    distance_norm_canberra<-c(distance_norm_canberra,dist2.norm)
  }
  return(list(distance_bray,distance_canberra,distance_norm_bray,distance_norm_canberra))
}



######################################################################################################################################
## Main program ##

#load arguments from command line
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}
array_id<-as.numeric(args[2])

# read in output from make_community.R
InComm<-list.files(pattern="comm.RData")

# loop through output, if multiple runs were performed
for (i in 1:length(InComm)) {
    load(InComm[i])
    model=strsplit(InComm[i],split="_")[[1]][2]
    basename=paste( substr(InComm[i],1,(nchar(InComm[i])-11)),"_it_",it,"_array_",array_id,sep="")
    ndat<-list()
    keyrun<-list()
    runA$Ab[is.finite(runA$Ab)==FALSE]<-0
    full_sites<-apply(runA$Ab,1,sum)!=0
    sites<-1:nrow(runA$Ab)
    availsites<-sites[full_sites]
    chosen_sites<-c()
    for(j in 1:it){
      site<-sample(availsites,1)
      chosen_sites<-c(chosen_sites,site)
      param<-list(r=runA$growth,a=runA$alpha,k=runA$carcap)
      keyruntmp<-LVMKeys(time=runA$time, param=param,init=runA$init[site,],AB=runA$Ab[site,])
      keyrun[[j]]<-keyruntmp
      keystoneDist<-keystone.distance(keyruntmp$Ab)
      #make networks
      g<-make_active_net(keyruntmp, nrow(keyruntmp$Ab), verb=F)
      g2<-graph.adjacency(keyruntmp$alpha,mode = "directed", weighted = TRUE)
      #do the stats
      cbtmp<-centrality.norm(g$graphs[[1]],type="betweenness",centralization=F)
      cctmp<-centrality.norm(g$graphs[[1]],type="closeness",centralization=F)
      cc_reltmp<-(centrality.norm(g$graphs[[1]],type="closeness",centralization=F))/max(centrality.norm(g$graphs[[1]],type="closeness",centralization=F))
      mdtmp<-(degree(g$graphs[[1]])-2)
      clusttmp<-transitivity(g$graphs[[1]],type="local")
      eigencentr<-eigen_centrality(g$graphs[[1]],directed=TRUE,weights=NA)$vector
      alphacentr<-alpha_centrality(g2)
      trans<-transitivity(g2,type = "weighted")
      ndattmp<-matrix(c(keystoneDist[[1]],keystoneDist[[2]],keystoneDist[[3]],keystoneDist[[4]],cbtmp,cctmp,cc_reltmp,mdtmp,clusttmp,eigencentr,alphacentr,trans),ncol=length(cbtmp),byrow=TRUE)
      rownames(ndattmp)<-c("bray.dist","canberra.dist","bray.dist.norm","canberra.dist.norm","cb","cc","cc_rel","md","clust","eigen.centr","alpha.centr","transitivity")
      ndat[[j]]<-ndattmp
      #list of lists
    }
runKey<-list(sites=chosen_sites,keyrun=keyrun,networkStats=ndat)
save(runKey, file=paste(basename,"_keystones.RData",sep=""))
}


############################################################################
### INPUT PARAMETERS

# it = number of iterations performed

## loads output of make_community.R (yourfilename_comm.RData); 
## this should therefore be in the same folder


############################################################################
### OUTPUT

# runKey

## $sites:

## $keyrun: (list of "it" entries, one per iteration; according to input parameter)
### $Ab = abundance table; columns are species, rows indicate which species was set to zero in this simulation run
### $init = initial abundances for each species and each simulated site (sampled from a uniform distribution, for details see documentation of seqtime::generateAbundances)
### $growth = growth rates per species (as specified by "mode" in input parameters)
### $carcap = carrying capacity (as specified by "mode" in input parameters)
### $alpha = interaction matrix (N x N; N = number of species simulated)

## $networkStats: (list of "it" entries, one per iteration; according to input parameter)
### each entry is a matrix containing network statistics for each species 
### (entries per species are identical in each matrix, except for the run where the species was set to 0)

### network statistics computed are:
#### bray.dist = Bray-Curtis dissimilarity (with vegdist); average BC-dissimilarity between communities without a particular species
#### canberra.dist = Canberra distance (with vegdist); average Canberra distance between communities without a particular species
#### bray.dist.norm = bray.dist normalized for the inital species abundances
#### canberra.dist.norm = canberra.dist normalized for the inital species abundances
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
# R CMD BATCH '--args it=1000 keystone_stat.R keystone_stat_log.txt
