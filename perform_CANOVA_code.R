library(quadprog)
library(ggplot2)
compute_statistic_twoside=function(data,index){
		index=factor(index)
		level=levels(index)
		Dmat <- matrix(0,length(level),length(level))
		diag(Dmat) <- table(index)
		Amat <- matrix(0,length(level)-1,length(level))
		for(i in 1:nrow(Amat)){
			Amat[i,c(i,i+1)]=c(1,-1)
			}
		Amat=t(Amat)
		bvec <- rep(0,length(level)-1)
		dvec=sum(data[index==level[1]])
		for(i in 2:length(level)){
			dvec=c(dvec,sum(data[index==level[i]]))
			}
		result=solve.QP(Dmat,dvec,Amat,bvec=bvec)$solution
		RSS1=sum((data[index==level[1]]-result[1])^2)
		for(i in 2:length(level)){
		RSS1=RSS1+sum((data[index==level[i]]-result[i])^2)
		}
		RSS0=sum((data[index==level[1]]-mean(data))^2)
		for(i in 2:length(level)){
		RSS0=RSS0+sum((data[index==level[i]]-mean(data))^2)
		}
		dist1=(RSS0-RSS1)/(RSS0/(length(data)-length(level)))

		Amat=-Amat
		result=solve.QP(Dmat,dvec,Amat,bvec=bvec)$solution
		RSS1=sum((data[index==level[1]]-result[1])^2)
		for(i in 2:length(level)){
		RSS1=RSS1+sum((data[index==level[i]]-result[i])^2)
		}
		RSS0=sum((data[index==level[1]]-mean(data))^2)
		for(i in 2:length(level)){
		RSS0=RSS0+sum((data[index==level[i]]-mean(data))^2)
		}
		dist2=(RSS0-RSS1)/(RSS0/(length(data)-length(level)))
		dists=max(dist1,dist2)
		return(dists)
	}


compute_statistic_oneside=function(data,index,side="increasing"){
		index=factor(index)
		level=levels(index)
		Dmat <- matrix(0,length(level),length(level))
		diag(Dmat) <- table(index)
		Amat <- matrix(0,length(level)-1,length(level))
		for(i in 1:nrow(Amat)){
			Amat[i,c(i,i+1)]=c(1,-1)
			}
		Amat=t(Amat)
		if(side=="decreasing"){
		Amat=-Amat
		}
		bvec <- rep(0,length(level)-1)
		dvec=sum(data[index==level[1]])
		for(i in 2:length(level)){
			dvec=c(dvec,sum(data[index==level[i]]))
			}
		result=solve.QP(Dmat,dvec,Amat,bvec=bvec)$solution
		RSS1=sum((data[index==level[1]]-result[1])^2)
		for(i in 2:length(level)){
		RSS1=RSS1+sum((data[index==level[i]]-result[i])^2)
		}
		RSS0=sum((data[index==level[1]]-mean(data))^2)
		for(i in 2:length(level)){
		RSS0=RSS0+sum((data[index==level[i]]-mean(data))^2)
		}
		dist1=(RSS0-RSS1)/(RSS0/(length(data)-length(level)))
		dists=dist1
		return(dists)
	}


dist_generate=function(numvec,iter=10000,side="both"){
result=c()
	for(itet1 in 1:iter){
		dat=rnorm(sum(numvec))
		index=rep(1,numvec[1])
		for(i in 2:length(numvec)){
		index=c(index,rep(i,numvec[i]))
		}
		index=factor(index)
		if(side=="both"){
		result=c(result,compute_statistic_twoside(dat,index))
		}else{
		result=c(result,compute_statistic_oneside(dat,index,side))
		}
	}
return(result)
}
test_with_one=function(dat,index,side="both",num_iter=10000,plotting=TRUE){
  dists=dist_generate(as.numeric(table(index)),num_iter,side=side)
  if(side=="both"){
  teststatistic=compute_statistic_twoside(dat,index)
  }else{
  teststatistic=compute_statistic_oneside(dat,index)
 }
 plot(hist(dists),xlim=c(0,max(c(dists,teststatistic))) ,main="distribution of null")
 abline(v=teststatistic,col="red")
 result=data.frame(numlevel=length(levels(index)),side=side,statistic=teststatistic,pvalue=(sum(teststatistic<=dists)+1)/(num_iter+1) )
 return(result)
}



test_bulk=function(dat,index,side="both",num_iter=10000){
	dists=dist_generate(as.numeric(table(index)),num_iter,side=side)
	stat=pval=c()
	if(side=="both"){
		for( i in 1:ncol(dat)){
		stat=c(stat,compute_statistic_twoside(dat[,i],index))
		}
	}else{
		for( i in 1:ncol(dat)){
		stat=c(stat,compute_statistic_oneside(dat[,i],index))
		}
	}

	for(i in 1:ncol(dat)){
		pval=c(pval,(sum(stat[i]<=dists)+1)/(num_iter+1))
	}	
	result=data.frame(idx=colnames(dat),side=side,statistic=stat,pvalue=pval,qvalue=p.adjust(pval,"fdr"))
	return(result)
}


test_with_one_cov=function(dat,index,covariate_mat,side="both",num_iter=10000,plotting=TRUE){
  covariate_mat=as.data.frame(covariate_mat)
  resi=residuals(lm(dat~.,data=covariate_mat))
  dists=dist_generate(as.numeric(table(index)),num_iter,side=side)
  if(side=="both"){
  teststatistic=compute_statistic_twoside(resi,index)
  }else{
  teststatistic=compute_statistic_oneside(resi,index)
 }
 plot(hist(dists),xlim=c(0,max(c(dists,teststatistic))) ,main="distribution of null")
 abline(v=teststatistic,col="red")
 result=data.frame(numlevel=length(levels(index)),side=side,statistic=teststatistic,pvalue=(sum(teststatistic<=dists)+1)/(num_iter+1) )
 return(result)
}



test_bulk_cov=function(dat,index,covariate_mat,side="both",num_iter=10000){
	dists=dist_generate(as.numeric(table(index)),num_iter,side=side)
	stat=pval=c()
	covariate_mat=as.data.frame(covariate_mat)
	if(side=="both"){
		for( i in 1:ncol(dat)){
	      resi=residuals(lm(dat[,i]~.,data=covariate_mat))
		stat=c(stat,compute_statistic_twoside(resi,index))
		}
	}else{
		for( i in 1:ncol(dat)){
	      resi=residuals(lm(dat[,i]~.,data=covariate_mat))
		stat=c(stat,compute_statistic_oneside(resi,index))
		}
	}

	for(i in 1:ncol(dat)){
		pval=c(pval,(sum(stat[i]<=dists)+1)/(num_iter+1))
	}	
	result=data.frame(idx=colnames(dat),side=side,statistic=stat,pvalue=pval,qvalue=p.adjust(pval,"fdr"))
	return(result)
}

perform_CANOVA=function(gene_dat,clinical_dat,side="both",num_iter=10000){
	gene_dat2=gene_dat[,c(1,match(colnames(gene_dat)[-1],clinical_dat$IndID,nomatch=-1)+1)]
	clinical_dat2=clinical_dat[match(clinical_dat$IndID,colnames(gene_dat2)[-1]),]
	dat=t(gene_dat2[,-1])
	colnames(dat)=gene_dat2[,1]
	index=clinical_dat2[,2]
	if(ncol(clinical_dat2)==2){
		if(ncol(dat)==1){
		return(test_with_one(dat,index,side=side,num_iter=num_iter))
		}else{
		return(test_bulk(dat,index,side=side,num_iter=num_iter))
		}
	}else{
	covariate_mat=clinical_dat2[,-c(1:2)]
		if(ncol(dat)==1){
		return(test_with_one_cov(dat,index,covariate_mat,side=side,num_iter=num_iter))
		}else{
		return(test_bulk_cov(dat,index,covariate_mat,side=side,num_iter=num_iter))
		}
	}
}

perform_CANOVA_external=function(gene_dat_path,clinical_dat_path,side="both",num_iter=10000){
	gene_dat=read.table(gene_dat_path,header=T,stringsAsFactor=F)
	clinical_dat=read.table(clinical_dat_path,header=T,stringsAsFactor=F)
	perform_CANOVA(gene_dat,clinical_dat,side=side,num_iter=num_iter)
	}


plot_fun=function(CANOVA_result,gene_dat,clinical_dat,top_num=1,cov_adjust=FALSE,path=""){
	CANOVA_result2=CANOVA_result[order(CANOVA_result$pvalue),]
	for(i in 1:top_num){
		png(paste(path,as.character(CANOVA_result2$idx[i]),"_top",top_num,".png",sep=""),width=600,height=600)
		dat2=gene_dat[which(gene_dat[,1]==as.character(CANOVA_result2$idx[i])),]
		gene_name=as.character(dat2[1,1])
		dat2=t(dat2[,-1])
		dat2=cbind(clinical_dat[,2],dat2)
		dat2=as.data.frame(dat2)
		colnames(dat2)=c("index","expr")
		dat2$index=factor(dat2$index)
		boxplot(dat2$expr~dat2$index,col=rainbow(length(unique(idx))) )
		dev.off()
	}
}


