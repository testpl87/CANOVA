library(quadprog)
compute_statistic_twoside=function(data,index){
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
