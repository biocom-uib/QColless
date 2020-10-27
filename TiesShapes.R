library(Zseq)
library(gmp)
library(ape)
library(igraph)
library(CollessLike)

## Computation of the ties

balance.indices2=function(tree){
  aux=c(round(CollessLike::colless.like.index(tree,f.size="ln",diss="MDM")/((1+log(exp(1)+2))/2)),
         round(CollessLike::colless.like.index(tree,f.size="ln",diss="var")/((1+log(exp(1)+2))^2/2)),
         CollessLike::sackin.index(tree),CollessLike::cophen.index(tree))
  return(aux)
}

collessQ = function(tree){
  return(CollessLike::colless.like.index(tree,f.size="ln",diss="var"))
}


are.tie = function(xx,yy) return(xx==yy)

simulations=5000
nf=18

  trees = list() 
  all.indices = list() 
  num.ties = c(0,0,0,0) 
  prob.ties = list() 
  all.num.ties = list()
  num.total.trees = list()
  for(n in 4:nf){
  dades = read.table(file=paste("IndicesTrees/indices_",n,".tsv",
                                 sep=""),header=TRUE)
    arbres = as.character(dades[,1])
    total.trees = length(arbres)
    total.pairs = total.trees*(total.trees-1)/2
    all.indices[[n]] = matrix(sapply(arbres,
                                   balance.indices2),ncol=4,byrow=T)
    num.ties[1]=sum(outer(all.indices[[n]][,1],
                        all.indices[[n]][,1],are.tie))
    num.ties[2]=sum(outer(all.indices[[n]][,2],
                        all.indices[[n]][,2],are.tie))
    num.ties[3]=sum(outer(all.indices[[n]][,3],
                        all.indices[[n]][,3],are.tie))
    num.ties[4]=sum(outer(all.indices[[n]][,4],
                          all.indices[[n]][,4],are.tie))
    num.ties = (num.ties-total.trees)/2
    prob.ties[[n]] = num.ties/total.pairs
    all.num.ties[[n]] = num.ties
    num.total.trees[[n]] = total.trees
    print(paste("Ties for n =",n," : ",
              paste(c("p_C=","p_C2","p_S=","p_Phi"),
                    round(prob.ties[[n]],4),collapse=", "),sep=""))
    print(paste("n=",n))
  }
  
  for (n in 19:20){
    dades = read.table(file=paste("IndicesTrees/indices_",n,".tsv",
                                  sep=""),header=TRUE)
    arbres.totals = as.character(dades[,1])
    arbres = sample(arbres.totals,simulations)
    total.trees = length(arbres)
    total.pairs = total.trees*(total.trees-1)/2
    all.indices[[n]] = matrix(sapply(arbres,
                                     balance.indices2),ncol=4,byrow=T)
    num.ties[1]=sum(outer(all.indices[[n]][,1],
                          all.indices[[n]][,1],are.tie))
    num.ties[2]=sum(outer(all.indices[[n]][,2],
                          all.indices[[n]][,2],are.tie))
    num.ties[3]=sum(outer(all.indices[[n]][,3],
                          all.indices[[n]][,3],are.tie))
    num.ties[4]=sum(outer(all.indices[[n]][,4],
                          all.indices[[n]][,4],are.tie))
    num.ties = (num.ties-total.trees)/2
    prob.ties[[n]] = num.ties/total.pairs
    all.num.ties[[n]] = num.ties
    num.total.trees[[n]] = total.trees
    print(paste("Ties for n =",n," : ",
                paste(c("p_C=","p_C2","p_S=","p_Phi"),
                      round(prob.ties[[n]],4),collapse=", "),sep=""))
    print(paste("n=",n))
  }
  


probabilities.ties = matrix(unlist(prob.ties),ncol=4,byrow=TRUE)
number.ties = matrix(unlist(all.num.ties),ncol=4,byrow=TRUE)
total.trees = matrix(unlist(num.total.trees),byrow=TRUE)

colnames(probabilities.ties)=c("Colless","Colless2","Sackin","Cophenetic")
colnames(number.ties)=c("Colless","Colless2","Sackin","Cophenetic")
rownames(probabilities.ties)=4:nf
rownames(number.ties)=4:nf
rownames(total.trees)=4:nf



## Plot of the ties

ties=read.table("TiesShapes.txt")
colnames(ties)=c("n","Colles","CollesQ","Sackin","Cophenetic","rQI","Total")
pr.colles=2*ties$Colles/(ties$Total*(ties$Total-1))
pr.collesq = 2*ties$CollesQ/(ties$Total*(ties$Total-1))
pr.sackin = 2*ties$Sackin/(ties$Total*(ties$Total-1))
pr.cofen = 2*ties$Cophenetic/(ties$Total*(ties$Total-1))
pr.qi = 2*ties$rQI/(ties$Total*(ties$Total-1))


table.ties = data.frame(4:20,pr.colles,pr.collesq,pr.sackin,pr.cofen,pr.qi)
names(table.ties)=c("n","Colless","Colless2","Sackin","Cophenetic","rQI")

x=4:20
require(ggplot2)

pdf(file="ties-shapes-paper.pdf")
ggplot(table.ties, aes(x)) +                   
  geom_line(aes(y=table.ties$Colless, colour="Colless")) +  
  geom_line(aes(y=table.ties$Colless2, colour="Q-Colless")) +
  geom_line(aes(y=table.ties$Sackin, colour="Sackin")) + 
  geom_line(aes(y=table.ties$Cophenetic, colour="Cophenetic")) +
  geom_line(aes(y=table.ties$rQI, colour="rQI")) +
  scale_y_continuous(trans='log',breaks=c(0,0.001,0.01,0.025,0.05,0.075),
                     labels=c("0","0.001","0.01","0.025","0.05","0.075")) +
  scale_colour_discrete("Method")

dev.off()
