source("perform_CANOVA_code.R")
genedat=matrix(rnorm(100*1000),nrow=1000,ncol=100)
genedat=data.frame(genedat)
genedat=cbind(paste("gene",1:1000,sep=""),genedat)
colnames(genedat)=c("geneID",paste("ind",1:100,sep=""))

idx=c(rep(1,25),rep(2,25),rep(3,25),rep(4,25))
covariate=data.frame(matrix(rnorm(100*5),nrow=100))
clinicaldat=cbind(paste("ind",1:100,sep=""),idx,covariate)
colnames(clinicaldat)=c("IndID","index",paste("cov",1:5,sep=""))

gene_dat=genedat
clinical_dat=clinicaldat
result=perform_CANOVA(genedat,clinicaldat,side="both",num_iter=1000)

data2=cbind(idx,as.numeric(gene_dat[1,-1]))
colnames(data2)=c("index","expr")
data2=data.frame(data2)
data2$index=factor(data2$index)

plot_fun(result,genedat,clinicaldat,10)


