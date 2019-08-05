##############################################################################################################################################
# SET OF FUNCTIONS TO COMPUTE eBIC and BIC criteria
# possible models : additive, additive+dominance, female+male, female+male+interaction
#
##note: require ASREML
#
##REQUIRED DATA & FORMAT
#
#PHENOTYPE - Y: a vector of length n, with names(Y)=individual names
#GENOTYPE  selec_XX: a list of length one, two or three matrices depending on the models
# matrices are n by mk matrix, where n=number of individuals, mk=number of selected SNPs,
#	with rownames()=individual names, and colnames()=SNP names
# additive: a single matrix, additive+dominance: two matrices
# female+male: two matrices with the female one first,
# female+male+interaction: three matrices with the female first, the male then the interaction
#KINSHIP  KK: a list of one, two or three matrices depending on the models
# additive: a n by n matrix, where n=number of individuals, with rownames()=colnames()=individual names
# additive+dominance: two n by n matrices, where n=number of individuals, with rownames()=colnames()=individual names
# female+male: a n.female by n.female matrix, with rownames()=colnames()=female names
#              a n.male by n.male matrix, with rownames()=colnames()=male names
# female+male+interaction: the same two matrices as the model female+male
#                 and a n by n matrix, where n=number of individuals, with rownames()=colnames()=individual names
#Factor - female: a factor of levels female names and length n, only for the last two models
#Factor - male: a factor of levels male names and length n, only for the last two models
#cofs:  a n by q matrix, where n=number of individuals, q=number of fixed effect,
# with rownames()=individual names and with column names, forbidden head of column names for this matrix "eff1_"
#
#Each of the previous data must be sorted in the same way, according to the individual name
#nb.tests: number of computed tests (total number of SNPs)
##FUNCTIONS USE
#save this file somewhere on your computer and source it!
#eBIC_allmodels(Y,selec_XX,KK,nb.tests,cofs=NULL,female=NULL,male=NULL)
# library(coxme)
# library(Matrix)
################################################################
# ind.covariances<-function(K,fact){ #already in mlmm script #it seems like this function is supposed to order K as the order of Y
#     RES<-matrix( unlist( sapply( fact, function(ii){
#         unlist( sapply( fact , function(jj){
#             il<-which(rownames(K)==levels(fact)[ii])
#             ic<-which(colnames(K)==levels(fact)[jj])
#             res<-K[il,ic]
#         }))
#     })),ncol=length(fact))
#     RES
# }
# #################################################################
# proj.matrix.sdp<-function(matrix){#already in mlmm script
# #projection of the matrix on definite positive matrix cone.
# cst<-0.000001
# mat.decomp<-eigen(matrix,symmetric=TRUE)
# valpp.mat<-mat.decomp$values
# valpp.mat[which(valpp.mat<cst)]<-cst  # transform too small or negative value
# valpp.proj<-valpp.mat
# res<-mat.decomp$vectors %*% (diag(valpp.proj)) %*% t(mat.decomp$vectors)
# colnames(res)<-colnames(matrix)
# rownames(res)<-rownames(matrix)
# res
# }
##################################################################
inv.matrix.sdp<-function(matrix){
#projection of the matrix on definite positive matrix cone and invertion.
	cst<-0.000001
	mat.decomp<-eigen(matrix,symmetric=TRUE)
	valpp.mat<-mat.decomp$values
	valpp.mat[which(valpp.mat<cst)]<-cst  # transform too small or negative value
	valpp.inv<-1/valpp.mat
	res<-mat.decomp$vectors %*% (diag(valpp.inv)) %*% t(mat.decomp$vectors)
	colnames(res)<-colnames(matrix)
	rownames(res)<-rownames(matrix)
	res
	}
###################################################################
my.model.lmekin<-function(cof.names,mod.random){
  res<-paste0("+",cof.names)
  res<-do.call(paste,c( as.list(res), sep=""))
  res<-as.formula( paste0 ("Y~0",res,mod.random) )
  res
}
###################################################################
my.random.lmekin<-function(eff.names){
  res<-paste0("+ (1|",eff.names,")")
  res<-do.call(paste,c( as.list(res), sep=""))
  res
}
##################################################################
new.BIC.lmekin <- function(lmekin.obj,lambda,nb.tests,fix,XX){ #following Gurka, TAS, 2006
  #06/2017 ML
  XX<-as.matrix(XX)
  logML<-lmekin.obj$loglik
  if(logML<0){
    fixsnp <- ncol(XX)
    nb.varcomp <-length(unlist(lmekin.obj$vcoef))+1
    rangXX<-Matrix::rankMatrix(XX)  #corrected to take the range of the fixed design matrix
    BIC <- -2*logML + log(nrow(XX)) * ( nb.varcomp + rangXX )
    ajout <-  2*lambda*lchoose( (nb.varcomp-1)*nb.tests, (fixsnp-fix) )
    eBIC <- BIC + ajout
    data.frame( BIC,ajout, eBIC, logML)
  }else{
    data.frame( NA,NA, NA, logML)
  }
}
########################################
lambda.calc<-function(n,nb.tests){
    if(n>=nb.tests){kk = 0}else{
      i = 0
      r = nb.tests
      while(r>n){
        i <- i+1
        r <- r/n
        kk=i
      }}
    if(kk == 0){
      lambda = 0
    }else{
      lambda = 1-(1/(2*kk))
    }
   lambda
}
###################################################################
#' @title Compute eBIC and BIC criteria
#' @param Y A numeric named vector where the names are individuals names and the values their phenotype.
#' @param selec_XX A list of length one, two or three matrices depending on the models. Use helper function \code{\link{frommlmm_toebic}} to get this argument.
#' @param KK a list of one, two or three matrices depending on the models
#'
#' - additive: a n by n matrix, where n=number of individuals, with rownames()=colnames()=individual names
#'
#' - additive+dominance: two n by n matrices, where n=number of individuals, with rownames()=colnames()=individual names
#'
#' - female+male: a n.female by n.female matrix, with rownames()=colnames()=female names and a n.male by n.male matrix, with rownames()=colnames()=male names
#'
#' - female+male+interaction: the same two matrices as the model female+male and a n by n matrix, where n=number of individuals, with rownames()=colnames()=individual names
#' @param nb.tests number of computed tests (total number of SNPs)
#' @param cofs  A n by q matrix, where n=number of individuals, q=number of fixed effect,
# with rownames()=individual names and with column names, forbidden head of column names for this matrix "eff1_"
#' @param female A factor of levels female names and length n, only for the last two models
#' @param male A factor of levels male names and length n, only for the last two models
#' @param lambda penalty used in the computation of the eBIC; if NULL, the default will be 1 - 1/(2k) with L=n^k where L=total number of SNPs (see function "lambda.calc")
#' @description Compute log likelihood, BIC and eBIC.
#'
#' The model with the smallest eBIC should be selected.
#' @return A matrix with a line for each mlmm step and 4 columns : BIC, ajout, eBIC_0.5 and LogL.
#' @template examples.donotrun
#' @export
eBIC_allmodels<-function(Y,selec_XX,KK,nb.tests,cofs=NULL,female=NULL,male=NULL,lambda=NULL) {
    #EDIT: adding special character support in marker names ( Olivier Guillaume 2018/07)
    for(ki in 1:length(selec_XX)) {
        stopifnot(anyDuplicated(colnames(selec_XX[[ki]])) == 0)
        if(ki == 1){
            newnames =  gsub("[^a-zA-Z0-9]","_",colnames(selec_XX[[ki]]))
            stopifnot(anyDuplicated(newnames) == 0)#different marker names become the same after special character being remplaced
            XX_mrk_names = structure( colnames(selec_XX[[ki]]), names = newnames)
        }else{
            stopifnot( colnames(selec_XX[[ki]]) == colnames(selec_XX[[1]]) )#checking if the markers are the same in the different matrices
        }
        colnames(selec_XX[[ki]]) = newnames
        #EDIT END
    }

    #edit by Olivier Guillaume 2018/02 centering and scaling Y because lmekins doesn't works well otherwise
    Ys = as.vector(scale(Y))
    names(Ys) = names(Y)
    Y = Ys
    rm(Ys)


    #edit by Olivier Guillaume 2018/02 sorting individuals alphabeticaly in Y, selec_XX, KK because of an lmekib bug
    indsSorted = sort(names(Y))
    i = match(indsSorted, names(Y))
    Y = Y[indsSorted]
    selec_XX = lapply(selec_XX, function(s_XX){
        m = as.matrix(s_XX[indsSorted,])
        colnames(m) = colnames(s_XX)
        m
    })
    if(is.null(female)){
        KK = lapply(KK, function(k){
            k[indsSorted,indsSorted]
        })
    }else{#edit by Olivier Guillaume 2018/04/04 this was forgotten
		female = female[i]
		male = male[i]
	}

nom.effet<-c("eff1","eff2","eff3")
n<-length(Y)
nb.effet<-length(selec_XX)
mk<-ncol(selec_XX[[1]])
for(ki in 1:nb.effet) {
	stopifnot(nrow(selec_XX[[ki]]) == n)
}
if(nb.effet > 1){
	for(ki in 2:nb.effet) {
		stopifnot(ncol(selec_XX[[ki]]) == mk)
	}
}
stopifnot(nrow(cofs) == n)
nb.level.byeffet<-rep(n,nb.effet)  #default for additive and dominance model
ind<-as.factor(names(Y))  #default for additive and dominance model
effet<-list()
for(ki in 1:nb.effet) {
	effet[[ki]]<-ind
}
stopifnot(length(KK) == nb.effet)
if(!is.null(female) & !is.null(male) ) {
	if(is.factor(female)==FALSE) female<-as.factor(female)
 	if(is.factor(male)==FALSE) male<-as.factor(male)
	n.female<-length(levels(female))
	n.male<-length(levels(male))
	if(nb.effet==2) nb.level.byeffet<-c(n.female,n.male)
	if(nb.effet==2) effet<-list(female,male)
	if(nb.effet==3) nb.level.byeffet<-c(n.female,n.male,n)
	if(nb.effet==3) effet<-list(female,male,ind)
}
for(ki in 1:nb.effet) {
	stopifnot(nrow(KK[[ki]]) == nb.level.byeffet[ki])
	stopifnot(ncol(KK[[ki]]) == nb.level.byeffet[ki])
}
for(ki in 1:nb.effet) {
	stopifnot(length(effet[[ki]]) == n)
}
#INTERCEPT
X0<-as.matrix(rep(1,n))
colnames(X0)<-"mu"
#kinship normalisation
KK.norm<-list()
for(ki in 1:nb.effet){
	n.temp<-nb.level.byeffet[ki]
	cst<-(n.temp-1)/sum((diag(n.temp)-matrix(1,n.temp,n.temp)/n.temp)*KK[[ki]] )
	KK.norm[[ki]]<-cst*KK[[ki]]
}
# rm(KK)
#contribution to individual covariances
KK.cov<-list()
for(ki in 1:nb.effet){
	KK.cov[[ki]]<-ind.covariances(KK.norm[[ki]],effet[[ki]])
	}
#projection on cone
KK.proj<-list()
for(ki in 1:nb.effet){
  KK.proj[[ki]]<-proj.matrix.sdp(KK.cov[[ki]])
}
names(KK.proj)<-nom.effet[1:nb.effet]

#build design matrix with correct order of columns
if(!is.null(cofs) ) X0<-cbind(X0,as.matrix(cofs))
snp.design<-NULL
for(ii in 1:ncol(selec_XX[[1]]) ){
	for(ki in 1:nb.effet){
		snp.design<-cbind(snp.design,selec_XX[[ki]][,ii])
	}
}
name.snp.col<-unlist(lapply(colnames(selec_XX[[1]]),function(xs){
	paste0(nom.effet[1:nb.effet],"_",xs)
	}))
design<-cbind(X0,snp.design)
colnames(design)<-c(colnames(X0),name.snp.col)

  if(is.null(lambda)){ ## edit by T. Flutre 16.07.19: possibility to pass lambda as arg
    ## lambda choice changed by Prune and Brigitte 23.02.16
    lambda<-lambda.calc(n,nb.tests)
  }
  #
  fix<-ncol(X0)
  snp<-nb.effet*ncol(selec_XX[[1]])
  #mod.random<-my.random.lmekin(names(effet))
	mod.random<-"+(1|ind)"
  result.eBIC <- t(matrix(unlist(lapply( seq(fix,(fix+snp),nb.effet), function(km){
    #
    xm<-my.model.lmekin(colnames(design)[1:km],mod.random)
    XX <- design[,1:km]
    #print(xm)
      tryCatch({
          #TODO remove suppressWarning once lmekin is updated. ie when depreciated rBind is replaced by rbind. (Olivier GUillaume 2018/05/04)
          suppressWarnings(res.lmekin<-coxme::lmekin(xm,data=data.frame(effet,Y, design,ind ), varlist=list(coxme::coxmeMlist(KK.proj) ) ,method="ML"))
          res.bic<- new.BIC.lmekin(res.lmekin,lambda,nb.tests,fix,XX)
          if(is.na(res.bic[1,1])){
            print(paste0("positive logML at step ",km))
          }
          return(res.bic)
      },
      error=function(cond){
        message(paste0("Problem at step ",km))
        message(cond)
        message(", Inf returned")
        return(rep(NA,4))
      })
})),nrow=4))
colnames(result.eBIC) <- c("BIC","ajout",paste("eBIC",lambda,sep="_"),"LogL")
fix.name<-do.call(paste, c(as.list(colnames(design)[1:fix],sep="+")))
rownames(result.eBIC) <- c(fix.name,colnames(selec_XX[[1]]))

#EDIT: adding special character support in marker names (Olivier Guillaume )
rownames(result.eBIC) = ifelse( rownames(result.eBIC) %in% names(XX_mrk_names), XX_mrk_names[rownames(result.eBIC)], rownames(result.eBIC))

#end of function
result.eBIC
}
