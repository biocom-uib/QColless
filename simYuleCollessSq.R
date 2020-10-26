library(mvSLOUCH)
library(apTreeshape)
library(ape)

flimdistDiscDirectRecursion.QCICI<-function(num.iter=10,popsize=10000,Y0=0,b_quad=TRUE){
    replicate(popsize,fDiscDirectDrawRecursion.QCICI(num.iter,Y0,b_quad))
}

fCtauDiscQCI<-function(x){1+6*x*x-6*x} ## for QCI
fCtauDiscCI<-function(x){x*log(x)+(1-x)*log(1-x)+1-2*min(x,1-x)} ## for Colless index

fDiscDirectDrawRecursion.QCICI<-function(num.iter=10,Y0=0,b_quad){
## Algorithm 3 from K. Bartoszek. Exact and approximate limit behaviour of the Yule tree's cophenetic index. Math. Biosci., 303:26–45, 2018.
    res<-0
    if (num.iter==0){
	Y1<-Y0
	Y2<-Y0
    }
    else{
	Y1<-fDiscDirectDrawRecursion.QCICI(num.iter-1,Y0,b_quad)
	Y2<-fDiscDirectDrawRecursion.QCICI(num.iter-1,Y0,b_quad)	
    }
    tau<-runif(1)
    if (b_quad){res<-((tau)^2)*Y1+((1-tau)^2)*Y2+fCtauDiscQCI(tau)}
    else{res<-((tau))*Y1+((1-tau))*Y2+fCtauDiscCI(tau)}
    res
}


f_getCollessSq<-function(phyltree){
## implementation of the definition of the Quadratic Colless index
    print(i_currenttree) ## counter which tree is being done at the moment
    i_currenttree<<-i_currenttree+1
    if (!inherits(phyltree,"phylo")){phyltree<-as.phylo(phyltree)}
    phyltree<-mvSLOUCH::phyltree_paths(phyltree) ## enhance the tree with extra fields in order to get number of "left" and "right" tip descendants
    lDesc<-as.list(rep(0,phyltree$Ntips-1))
    for (i_tip in 1:phyltree$Ntips){
	for (i_intnode in (rev(phyltree$path.from.root[[i_tip]]$nodes)[-1])){## the final tip node is also in this vector
	    lDesc[[i_intnode-phyltree$Ntips]]<-lDesc[[i_intnode-phyltree$Ntips]]+1
	}	
    }
    iCollessSq<-0
    for (i_intnode in phyltree$internal_nodes_index){
	vLRchildren<-phyltree$edge[which(phyltree$edge[,1]==i_intnode),2]
	Ldesc<-1
	Rdesc<-1
	if (length(vLRchildren)==2){
	    if (vLRchildren[1]>(phyltree$Ntips)){Ldesc<-lDesc[[vLRchildren[1]-phyltree$Ntips]]}
	    if (vLRchildren[2]>(phyltree$Ntips)){Rdesc<-lDesc[[vLRchildren[2]-phyltree$Ntips]]}
	}
	iCollessSq<-iCollessSq+(Ldesc-Rdesc)^2
    }
    iCollessSq
}

b_replicate_study<-TRUE  ## set to FALSE if an independent run is wanted and not a replication
s_randomseedfile<-"simYuleCollessSq_randomseed.RData"
s_outputfile<-"simYuleCollessSq_simulatedvalues.RData"
if (b_replicate_study && file.exists("s_randomseedfile")){
    load("simYuleCollessSq.RData")
    RNGkind(kind = RNG_kind[1], normal.kind = RNG_kind[2], sample.kind = RNG_kind[3])
    RNGversion(RNG_version)  
}else{
    rexp(1)
    RNG_version<-getRversion()
    RNG_kind<-RNGkind()
    save(.Random.seed,RNG_kind,RNGversion,file=s_randomseedfile)
}

numtrees<-10000


##=====================================================================

num.iter<-10 ## number of iterations for the recursive algorithm to simulate from the finite step approximation to the limit distribution
Y0 <- 0 ## initial value of the recursion QCI, CI for tree with two tips

## simulate from finite step approximation to limit distribution of the Quadratic Colless index
limapproxdens_recursion_10iter_QCI<-flimdistDiscDirectRecursion.QCICI(num.iter=num.iter,popsize=numtrees,b_quad=TRUE)
pdf("YuleQCItheor_10iter.pdf");hist(limapproxdens_recursion_10iter_QCI,breaks=100,col="black",freq=FALSE,xlab="Yule QCi",ylab="",cex.lab=1.5,cex.axis=1.5,main="",ylim=c(0,1),xlim=c(-2.5,3.5));dev.off()


## simulate from finite step approximation to limit distribution of the Colless index
limapproxdens_recursion_10iter_CI<-flimdistDiscDirectRecursion.QCICI(num.iter=num.iter,popsize=numtrees,b_quad=FALSE)
pdf("YuleCItheor_10iter.pdf");hist(limapproxdens_recursion_10iter_CI,breaks=100,col="black",freq=FALSE,xlab="Yule Colless index",ylab="",cex.lab=1.5,cex.axis=1.5,main="",ylim=c(0,1),xlim=c(-2.5,3.5));dev.off()


##=====================================================================

## histogram of the normalized QC index for Yule trees
ntip<-1000 
i_currenttree<-1
vCollessSqYule<-sapply(apTreeshape::rtreeshape(numtrees,tip.number=ntip,model="yule"),FUN=f_getCollessSq,simplify=TRUE)
meanCollessSqYule<-ntip*(ntip+1)-2*ntip*sum(1/(ntip:1))
sdleadingtermCollessSqYule<-ntip^2
vCollessSqYulenorm<-(vCollessSqYule-meanCollessSqYule)/sdleadingtermCollessSqYule
pdf("YuleQCIsim.pdf");hist(vCollessSqYulenorm,breaks=100,col="black",freq=FALSE,xlab="Yule squared Colless index",ylab="",cex.lab=1.5,cex.axis=1.5,main="",ylim=c(0,1),xlim=c(-2.5,3.5));dev.off()

save(limapproxdens_recursion_10iter_QCI,limapproxdens_recursion_10iter_CI,vCollessSqYule,vCollessSqYulenorm,file=s_outputfile)
