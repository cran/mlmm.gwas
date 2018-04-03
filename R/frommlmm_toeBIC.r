#' @title Helper that create the selec_XX argument of eBIC_allmodels()
#' @description Helper function that create the selec_XX argument of \code{\link{eBIC_allmodels}} from the output of \code{\link{mlmm_allmodels}}.
#' @param XX A list of length one, two or three matrices depending on the model. Matrices are n by m matrix, where n=number of individuals, m=number of SNPs, with rownames(X)=individual names, and colnames(X)=SNP names.
#' Use the same XX you used with the \code{\link{mlmm_allmodels}} function
#' @param res.mlmm Output from the \code{\link{mlmm_allmodels}} function
#' @seealso \code{\link{mlmm_allmodels}} \code{\link{eBIC_allmodels}}
#' @template examples.donotrun
#' @export
frommlmm_toebic<-function(XX, res.mlmm){
    n.step<-length(res.mlmm)
    id<-grep("selec_",names(res.mlmm[[n.step]]) )
    names.snp<-names(res.mlmm[[n.step]])[id]
    names.snp<-unlist(sapply(names.snp,function(x){
    	unlist(strsplit(x,"selec_"))[2]
    	}))
    #add the last detected snp
    if(n.step>2) {
    	id<-which( res.mlmm[[n.step]][-(1:(n.step-2))]==min( res.mlmm[[n.step]][-(1:(n.step-2))] ,na.rm=TRUE) )
    	last.snp<-names( res.mlmm[[n.step]])[-(1:(n.step-2))][id]
    	} else {
    	id<-which(res.mlmm[[n.step]]==min( res.mlmm[[n.step]],na.rm=TRUE ) )
    	last.snp<-names( res.mlmm[[n.step]])[id]
    	}
    names.snp<-c(names.snp,last.snp)
    nb.effet<-length(XX)
    selec_XX<-list()
    for(ki in 1:nb.effet){
    	selec_XX[[ki]]<-matrix(unlist(lapply(names.snp,function(xj){
    	XX[[ki]][,which(colnames(XX[[ki]])==xj)]
    	})),ncol=length(names.snp),byrow=FALSE)
    	}
    for(ki in 1:nb.effet){
    	colnames(selec_XX[[ki]])<-names.snp
    	rownames(selec_XX[[ki]])<-rownames(XX[[1]])
    	}
    #end of function
    res<-selec_XX
}

