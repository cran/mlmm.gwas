###adapted from MLMM - Multi-Locus Mixed Model
#' @import stats
#' @import grDevices

##############################################################################################################################################
# SET OF FUNCTIONS TO CARRY GWAS CORRECTING FOR POPULATION STRUCTURE WHILE INCLUDING
# COFACTORS THROUGH A FORWARD REGRESSION APPROACH
# possible models : additive, additive+dominance, female+male, female+male+interaction
#
##REQUIRED DATA & FORMAT
#
#PHENOTYPE - Y: a vector of length n, with names(Y)=individual names
#GENOTYPE  XX: a list of length one, two or three matrices depending on the models
# matrices are n by m matrix, where n=number of individuals, m=number of SNPs,
#	with rownames(X)=individual names, and colnames(X)=SNP names
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
# with rownames()=individual names and with column names, forbidden head of column names for this matrix "eff1_" and usage of
# special characters as "*","/","&"
#
#Each of the previous data must be sorted in the same way, according to the individual name
##FUNCTIONS USE
#save this file somewhere on your computer and source it!
#
#mlmm_allmodels(Y,XX,KK,nbchunks,maxsteps,cofs=NULL,female=NULL,male=NULL)
#
#nbchunks: an integer defining the number of chunks of matrices to run the analysis, allows to decrease the memory usage ==> minimum=2, increase it if you do not have enough memory
#maxsteps: maximum number of steps desired in the forward approach. The forward approach breaks automatically once the pseudo-heritability is close to 0,
#however to avoid doing too many steps in case the pseudo-heritability does not reach a value close to 0, this parameter is also used.
#It's value must be specified as an integer >= 3
# library(Matrix)
# library(sommer)
##################################################################
Prep.forRSS<-function(istep,InvV,Y,cof_fwd,n){
    #this function computes matrices to accelerate the test computations
    Y_t<-crossprod( InvV ,Y)
    cof_fwd_t<-crossprod( InvV,cof_fwd[[(istep-1)]] )
    Res_H0<-summary(lm(Y_t~0+cof_fwd_t))$residuals
    Q_<-qr.Q(qr(cof_fwd_t))
    M<-InvV %*% (diag(n)-tcrossprod(Q_,Q_))
    res<-list(M=M,cof_fwd_t=cof_fwd_t,Res_H0=Res_H0)
    rm(Y_t,Q_)
    res
}
##################################################################
find.names<-function(xx){
    res<-unique(unlist( strsplit ( xx,"eff1_") ))
    res
}
##################################################################
RSS.forSNP<-function(istep,list.Prep,nbchunks,XX){
    #this function computes the residual sum of squares
    RSS<-list()
    XXt<-list()
    names.cof_fwd<-find.names( colnames(list.Prep$cof_fwd_t) )
    coltokeep<-which(!colnames(XX[[1]]) %in% names.cof_fwd )
    if(round(length(coltokeep)/nbchunks) !=0 ) {
        for (j in 1:(nbchunks-1)) {
            coltoget<-((j-1)*round(length(coltokeep)/nbchunks)+1):(j*round(length(coltokeep)/nbchunks))
            for(ki in 1:length(XX) ){
                XXt[[ki]]<-crossprod( list.Prep$M ,as.matrix(XX[[ki]][,coltokeep])[,coltoget] )
            }
            RSS[[j]]<-unlist(sapply(1:ncol(XXt[[1]]),function(iij){
                x<-NULL
                for(ki in 1:length(XX) ){
                    x<-cbind(x,XXt[[ki]][,iij])
                }
                temp<-lsfit(x,list.Prep$Res_H0,intercept = FALSE)
                res<-c( sum(temp$residuals^2),temp$qr$rank)
                res
            }
            ))
        }
    }
    if( ( (nbchunks-1)*round(length(coltokeep)/nbchunks)+1) <= length(coltokeep) ) {
        coltoget<-( ( (nbchunks-1)*round( length(coltokeep)/nbchunks)+1):length(coltokeep))
        for(ki in 1:length(XX) ){
            XXt[[ki]]<-crossprod( list.Prep$M ,as.matrix(XX[[ki]][,coltokeep])[,coltoget] )
        }
        RSS[[nbchunks]]<-unlist(sapply(1:ncol(XXt[[1]]),function(iij){
            x<-NULL
            for(ki in 1:length(XX) ){
                x<-cbind(x,XXt[[ki]][,iij])
            }
            temp<-lsfit(x,list.Prep$Res_H0,intercept = FALSE)
            res<-c( sum(temp$residuals^2), temp$qr$rank)
            res
        }
        ))
    }
    rm(XXt)
    RES<-matrix(unlist(RSS),ncol=2,byrow=T)
    rownames(RES)<-colnames(XX[[1]])[coltokeep]
    colnames(RES)<-c("RSS","rank")
    RES
}
####################################################################
my.pval.aov<-function(addcof_fwd,cof_fix,InvV,Y,XX){
    #this function compute p value of fixed effect
    Y_t<-crossprod(InvV,Y)
    reg<-list()
    reg<-lapply(2:length(addcof_fwd),function(ii){
        res<-NULL
        for(ki in 1:length(XX) ) {
            res<-cbind(res,XX[[ki]][,colnames(XX[[1]]) %in% addcof_fwd[[ii]] ] )
        }
        res<-crossprod(InvV,res)
        res
    })
    fix<-crossprod(InvV,cof_fix)
    model<-"Y_t~0+fix"
    for(i in 1:(length(addcof_fwd)-1) ){
        model<-paste0(model,"+reg[[",i,"]]")
    }
    model<-as.formula(model)
    temp<-drop1(aov(model),test="F")  #type II Sum of squares
    rownames.wob<-unlist(sapply(rownames(temp),function(xn){unlist(strsplit(xn," "))}))
    if("fix" %in% rownames.wob){res<-temp$"Pr(>F)"[-c(1,2)]} else { res<-temp$"Pr(>F)"[-1]}
    names(res)<-paste0("selec_",unlist(addcof_fwd)[2:length(addcof_fwd)])
    res
}
################################################################
ind.covariances<-function(K,fact){
    RES<-matrix( unlist( sapply( fact, function(ii){
        unlist( sapply( fact , function(jj){
            il<-which(rownames(K)==levels(fact)[ii])
            ic<-which(colnames(K)==levels(fact)[jj])
            res<-K[il,ic]
        }))
    })),ncol=length(fact))
    RES
}
##################################################################
proj.matrix.sdp<-function(matrix){
    #projection of the matrix on definite positive matrix cone.
    cst<-0.000001
    mat.decomp<-eigen(matrix,symmetric=TRUE)
    valpp.mat<-mat.decomp$values
    valpp.mat[which(valpp.mat<cst)]<-cst  # transform too small or negative value
    valpp.proj<-valpp.mat
    res<-mat.decomp$vectors %*% (diag(valpp.proj)) %*% t(mat.decomp$vectors)
    colnames(res)<-colnames(matrix)
    rownames(res)<-rownames(matrix)
    res
}
###################################################################
mixedModelWithSommer = function(Y, KK.proj, effet, cof=NULL){
    data = data.frame(id=names(Y), Y=as.vector(Y), effet, cof)
    random = formula(paste("~",paste0("g(",names(KK.proj),")",collapse="+")))
    fixed = formula(sub("mu", "1", paste("Y~",paste(colnames(cof), collapse = "+"))))
    currentDevice = dev.cur()
    res = sommer::mmer2(fixed, random, data=data, G=KK.proj, silent=TRUE)#NOTE: opens an unwanted empty graphical window (probably a bug of sommer).
    if(dev.cur() != 1 && dev.cur() != currentDevice) dev.off()#EDIT: close the unwanted window. temporary fix, it would be better to not opens this window in the first time (OG 2018/07/17)
    res
}
###################################################################
#' @title Multi-Locus Mixed-Model
#'
#' @description
#' Carry GWAS correcting for population structure while including cofactors through a forward regression approach.
#'
#' possible models: additive, additive+dominance, female+male, female+male+interaction
#'
#' For additive model, look at the example below or at this \href{../doc/gwas-manual.html}{vignette}. For other models, read \href{https://link.springer.com/article/10.1007/s00122-017-3003-4}{Bonnafous et al. (2017)}.
#'
#' @param Y A numeric named vector where the names are individuals' names and the values their phenotype. The names of Y will be matched to the row names of X.
#' @param XX A list of length one, two or three matrices depending on the model. Matrices are n by m matrix, where n=number of individuals, m=number of SNPs, with rownames(X)=individual names, and colnames(X)=SNP names.
#'
#' - additive: a single matrix
#'
#' - additive+dominance: two matrices
#'
#' - female+male: two matrices with the female one first
#'
#' - female+male+interaction: three matrices with the female first, the male then the interaction
#' @param KK a list of one, two or three matrices depending on the models
#'
#' - additive: a n by n matrix, where n=number of individuals, with rownames()=colnames()=individual names
#'
#' - additive+dominance: two n by n matrices, where n=number of individuals, with rownames()=colnames()=individual names
#'
#' - female+male: a n.female by n.female matrix, with rownames()=colnames()=female names and a n.male by n.male matrix, with rownames()=colnames()=male names
#'
#' - female+male+interaction: the same two matrices as the model female+male and a n by n matrix, where n=number of individuals, with rownames()=colnames()=individual names
#' @param nbchunks An integer defining the number of chunks of matrices to run the analysis, allows to decrease the memory usage. minimum=2, increase it if you do not have enough memory
#' @param maxsteps An integer >= 3. Maximum number of steps desired in the forward approach. The forward approach breaks automatically once the pseudo-heritability is close to 0, however to avoid doing too many steps in case the pseudo-heritability does not reach a value close to 0, this parameter is also used.
#' @param cofs A n by q matrix, where n=number of individuals, q=number of fixed effect, with rownames()=individual names and with column names, forbidden head of column names for this matrix "eff1_" and usage of special characters as "*","/","&"
#' @param female A factor of levels female names and length n, only for the last two models
#' @param male A factor of levels male names and length n, only for the last two models
#' @details Each of the data arguments must be sorted in the same way, according to the individual name.
#' @return a list with one element per step of the forward approach. Each element of this list is a named vector of p-values, the names are the names of the markers, with "selec_" as prefix for the markers used as fixed effects.
#' @seealso \code{\link{manhattan.plot}}
#' @template examples
#' @export
mlmm_allmodels<-function(Y, XX, KK, nbchunks=2, maxsteps=20, cofs=NULL, female=NULL, male=NULL) {
    nom.effet<-c("eff1","eff2","eff3")
    n<-length(Y)
    nb.effet<-length(XX)
    m<-ncol(XX[[1]])
    for(ki in 1:nb.effet) {
        stopifnot(nrow(XX[[ki]]) == n)
        #EDIT OG 2018/07/10: adding support of marker names in the form "chr3:10179093_T/G" (incompatible when used directly in formulas)
        stopifnot(anyDuplicated(colnames(XX[[ki]])) == 0)
        if(ki == 1){
            newnames =  gsub("[^a-zA-Z0-9]","_",colnames(XX[[ki]]))
            stopifnot(anyDuplicated(newnames) == 0)#different marker names become the same after special character being remplaced
            XX_mrk_names = structure( colnames(XX[[ki]]), names = newnames)
        }else{
            stopifnot( colnames(XX[[ki]]) != colnames(XX[[1]]) )#checking if the markers are the same in the different matrices
        }
        colnames(XX[[ki]]) = newnames
        #EDIT END
    }
    if(nb.effet > 1){
        for(ki in 2:nb.effet) {
            stopifnot(ncol(XX[[ki]]) == m)
        }
    }
    stopifnot(nrow(cofs) == n)
    nb.level.byeffect<-rep(n,nb.effet)  #default for additive and dominance model
    ind<-as.factor(names(Y)) #default for additive and dominance model
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
        if(nb.effet==2) nb.level.byeffect<-c(n.female,n.male)
        if(nb.effet==2) effet<-list(female,male)
        if(nb.effet==3) nb.level.byeffect<-c(n.female,n.male,n)
        if(nb.effet==3) effet<-list(female,male,ind)
    }
    for(ki in 1:nb.effet) {
        stopifnot(nrow(KK[[ki]]) == nb.level.byeffect[ki])
        stopifnot(ncol(KK[[ki]]) == nb.level.byeffect[ki])
    }
    for(ki in 1:nb.effet) {
        stopifnot(length(effet[[ki]]) == n)
    }
    stopifnot(nbchunks >= 2)
    stopifnot(maxsteps >= 3)
    #INTERCEPT
    X0<-as.matrix(rep(1,n))
    colnames(X0)<-"mu"
    #kinship normalisation
    KK.norm<-list()


    for(ki in 1:nb.effet){
        n.temp<-nb.level.byeffect[ki]
        cst<-(n.temp-1)/sum((diag(n.temp)-matrix(1,n.temp,n.temp)/n.temp)*KK[[ki]] )
        KK.norm[[ki]]<-cst*KK[[ki]]
    }

    #contribution to individual covariances
    KK.cov<-list()
    for(ki in 1:nb.effet){
        KK.cov[[ki]]<-ind.covariances(KK.norm[[ki]],effet[[ki]])
    }
    #change for lmekin, all random effects must be ind
    for(ki in 1:nb.effet) {
        effet[[ki]]<-ind
    }
    names(effet)<-nom.effet[1:nb.effet]
    #projection on cone
    KK.proj<-list()
    for(ki in 1:nb.effet){
        KK.proj[[ki]]<-proj.matrix.sdp(KK.cov[[ki]])

        colnames(KK.proj[[ki]]) = names(Y)
        rownames(KK.proj[[ki]]) = names(Y)
    }
    names(KK.proj)<-nom.effet[1:nb.effet]
    ############################################to discard
    #step 0 : NULL MODEL
    addcof_fwd<-list()
    addcof_fwd[[1]]<-'NA'
    #
    cof_fwd<-list()
    if(is.null(cofs) ) {cof_fwd[[1]]<-X0 } else{ cof_fwd[[1]]<-cbind(X0,as.matrix(cofs) )}
    mod_fwd<-list()
    mod_fwd[[1]] = mixedModelWithSommer(Y, KK.proj, effet, cof_fwd[[1]])

    herit_fwd <- list()
    herit_fwd[[1]] = 1-mod_fwd[[1]]$sigma['units']/sum(mod_fwd[[1]]$sigma)

    V = mod_fwd[[1]]$sigma['units']*diag(n)


    for(ki in 1:nb.effet){
        V = V + mod_fwd[[1]]$sigma[[ki]]*KK.cov[[ki]]
    }
    InvV <- solve(chol(V))
    df1 <- list()
    df1[[1]] <- Matrix::rankMatrix(cof_fwd[[1]])
    pval<-list()
    pval[[1]]<-NA
    fwd_lm<-list()
    pval_cof_fwd<-list()
    pval_cof_fwd[[1]]<-NULL
    pval_cof_fwd[[2]]<-NULL
    cat('null model done! pseudo-h=',round(herit_fwd[[1]],3),'\n')

    #
    for (i in 2:(maxsteps)) {
        if (herit_fwd[[i-1]] < 0.01) break else {

            list.Prep<-Prep.forRSS(i,InvV,Y,cof_fwd,n)
            RSS_H1.rank<-RSS.forSNP(i,list.Prep,nbchunks,XX)
            RSS_H1<-RSS_H1.rank[,1]
            RSS_H0<-sum( list.Prep$Res_H0^2 )
            df2<-n-RSS_H1.rank[,2]-df1[[i-1]]
            df1.test<-RSS_H1.rank[,2]
            Ftest<-(rep( RSS_H0,length( RSS_H1 ))/  RSS_H1-1)*df2/df1.test
            Ftest[is.na(Ftest)]<-0
            pval[[i]]<-pf(Ftest,df1.test,df2,lower.tail=FALSE)

            addcof_fwd[[i]]<-names(which(RSS_H1==min(RSS_H1))[1])
            cof_fwd[[i]]<-cof_fwd[[ (i-1) ]]
            for(ki in 1:nb.effet){
                cof_fwd[[i]]<-cbind(cof_fwd[[i]], XX[[ki]][,colnames(XX[[1]]) %in% addcof_fwd[[i]] ])
            }
            colnames(cof_fwd[[i]])[( ncol(cof_fwd[[i]])-(nb.effet-1) ):ncol(cof_fwd[[i]]) ]<-c(paste0( nom.effet[1:nb.effet],"_",addcof_fwd[[i]]) )

            mod_fwd[[i]] = mixedModelWithSommer(Y, KK.proj, effet, cof_fwd[[i]])

            df1[[i]] = Matrix::rankMatrix(cof_fwd[[i]])
            herit_fwd[[i]] = 1-mod_fwd[[i]]$sigma['units']/sum(mod_fwd[[i]]$sigma)

            rm(list.Prep)
            cat('step ',i-1,' done! pseudo-h=',round(herit_fwd[[i]],3),'model: Y ~',paste0(sub("^eff1_","",colnames(cof_fwd[[i]])),collapse=" + "),'\n')
        }
        ##get pval at each forward step for selected SNP

        V = mod_fwd[[i]]$sigma['units']*diag(n)
        for(ki in 1:nb.effet){
            V = V + mod_fwd[[i]]$sigma[[ki]]*KK.cov[[ki]]
        }
        InvV <- solve(chol(V))
        pval_cof_fwd[[i+1]]<-my.pval.aov(addcof_fwd,cof_fwd[[1]],InvV,Y,XX)
        pval[[i]]<-c(pval_cof_fwd[[i]],pval[[i]])
    }
    #end of function

    #EDIT putting back the markers names with special caracters
    if(length(pval)>=2){
        for(i in 2:length(pval)){
            isSelec = grepl('^selec_',names(pval[[i]]))
            newnames = XX_mrk_names[ sub("^selec_","",names(pval[[i]])) ]
            newnames = ifelse(isSelec, paste0("selec_",newnames), newnames)
            names(pval[[i]]) = newnames
        }
    }


    pval
}


