#' @title Helper function that create the selec_XXclass argument of the Estimation_allmodels function
#' @description Function that create the selec_XXclass argument of the \code{\link{Estimation_allmodels}} function from the output of the \code{\link{eBIC_allmodels}} function.
#' @param XX A list of length one, two or three matrices depending on the model. Matrices are n by m matrix, where n=number of individuals, m=number of SNPs, with rownames(X)=individual names, and colnames(X)=SNP names.
#'
#' - additive: a single matrix
#'
#' - additive+dominance: two matrices
#'
#' - female+male: two matrices with the female one first
#'
#' - female+male+interaction: three matrices with the female first, the male then the interaction
#' @param res.eBIC output of the \code{\link{eBIC_allmodels}} function
#' @seealso \code{\link{eBIC_allmodels}} \code{\link{Estimation_allmodels}}
#' @template examples.donotrun
#' @export
fromeBICtoEstimation<-function(XX,res.eBIC){
##########################################
#GENOTYPE  XX: a list of length one or two matrices depending on the models
# matrices are n by m matrix, where n=number of individuals, m=number of SNPs,
#	with rownames(X)=individual names, and colnames(X)=SNP names
# To get a three class genotypes (00,01+10,11), only a single matrix is necessary that of the additive model
# To get a four class haplotypes (00,10,01,11), two matrices those of the female+male model with the female one first
# this code assumes that the genotype matrices have been centered but works if they are not
#
#res.eBIC: output of eBIC function, 4 column matrix
#
#output is NULL if no SNP is selected
######################################
nb.effet<-length(XX)
eBIC<-res.eBIC[,3]
il<-which( eBIC==min(eBIC, na.rm = TRUE ))# ,na.rm added by lgody on 2018/01/23
snp.selec<-NULL
selec_XXclass<-NULL
if(il>1) {selec_XX<-list()
	snp.selec<-rownames(res.eBIC)[2:il]
	for(ii in 1:nb.effet){
		#EDIT by Olivier Guillaume (2018/07/17): conserve order of eBIC table
	    #selec_XX[[ii]]<-XX[[ii]][,colnames(XX[[ii]])%in%snp.selec]
	    selec_XX[[ii]]<-XX[[ii]][,match(snp.selec, colnames(XX[[ii]]))]

		selec_XX[[ii]] = as.matrix(selec_XX[[ii]]) ## add by Prune 15.06.16
		}
	for(ii in 1:nb.effet){ stopifnot( ncol(selec_XX[[ii]]) == length(snp.selec) ) }
	if(nb.effet == 1){
		selec_XX<-as.matrix(selec_XX[[1]])
		selec_XXclass<-data.frame(lapply(1:length(snp.selec) ,function(ic){
			mini<-min( selec_XX[,ic] )
			maxi<-max( selec_XX[,ic] )
			trans<-unlist(lapply(selec_XX[,ic] , function(xx){
				res<-"01|10"
				if(xx == mini) res<-"00"
				if(xx == maxi) res<-"11"
				res
				}))  #the original codage is (0,1,2), 0=XRQ
			trans<-as.factor(trans)
			}))
		colnames(selec_XXclass)<-snp.selec
	}
	if(nb.effet == 2){
		minusfreq<-matrix(NA,ncol=length(snp.selec),nrow=nb.effet)
		for(ii in 1:nb.effet){
			minusfreq[ii,]<-apply(selec_XX[[ii]],2,min)
			}  #the original codage is (0,1), 0=XRQ
		selec_XX_t<-list()
		for(ii in 1:nb.effet){
		selec_XX_t[[ii]]<-selec_XX[[ii]]-matrix(rep( minusfreq[ii,],nrow(selec_XX[[1]]) ),ncol=length(snp.selec),nrow=nrow(selec_XX[[1]]),byrow=TRUE)
		}
		selec_XXclass<-data.frame(lapply(1:length(snp.selec) ,function(ic){
			trans<-unlist(lapply(1:nrow(selec_XX[[1]]),function(il){
				paste0( selec_XX_t[[1]][il,ic],selec_XX_t[[2]][il,ic] )
				}))
			trans<-as.factor(trans)
		}))
		colnames(selec_XXclass)<-snp.selec
	}
}
#end of function
selec_XXclass
}
