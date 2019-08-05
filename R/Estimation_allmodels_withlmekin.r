##############################################################################################################################################
# SET OF FUNCTIONS TO COMPUTE etimated effects
#
#note: require ASREML, multcomp and multcompView
#
##REQUIRED DATA & FORMAT
#
#PHENOTYPE - Y: a vector of length n, with names(Y)=individual names
#GENOTYPE  selec_XXclass: a n by mk data.frame of factors
# with rownames()=individual names, and colnames()=mk selected SNP names
# additive+dominance: three levels factor
# female+male+interaction: four levels factor
#KINSHIP  KK: a list  two or three matrices depending on the models
# additive+dominance: two n by n matrices, where n=number of individuals, with rownames()=colnames()=individual names
# female+male+interaction: a n.female by n.female matrix, with rownames()=colnames()=female names
#              	           a n.male by n.male matrix, with rownames()=colnames()=male names
#                 and a n by n matrix, where n=number of individuals, with rownames()=colnames()=individual names
#Factor - female: a factor of levels female names and length n, only for the last two models
#Factor - male: a factor of levels male names and length n, only for the last two models
#cofs:  a n by q matrix, where n=number of individuals, q=number of fixed effet,
# with rownames()=individual names and with column names, forbidden head of column names for this matrix "eff1_"
#
#Each of the previous data must be sorted in the same way, according to the individual name
##FUNCTIONS USE
#save this file somewhere on your computer and source it!
#Estimation_allmodels(Y,selec_XXclass,KK,cofs=NULL,female=NULL,male=NULL)
#
# library(multcomp)
# library(multcompView)
# library(coxme)
# library(Matrix)
################################################################
# ind.covariances<-function(K,fact){#already in mlmm script
	# RES<-matrix( unlist( sapply( fact, function(ii){
		# unlist( sapply( fact , function(jj){
			# il<-which(rownames(K)==levels(fact)[ii])
			# ic<-which(colnames(K)==levels(fact)[jj])
			# res<-K[il,ic]
			# }))
		# })),ncol=length(fact))
	# RES
# }
##################################################################
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
###################################################################
# my.model.lmekin<-function(cof.names,mod.random){#already in eBic script.
  # res<-paste0("+",cof.names)
  # res<-do.call(paste,c( as.list(res), sep=""))
  # res<-as.formula( paste0 ("Y~0",res,mod.random) ) # 1 modified by Prune 27.07.17 reverted to 0 Gody + Mangin 2018.01.16
  # res
# }
###################################################################
# my.random.lmekin<-function(eff.names){
  # res<-paste0("+ (1|",eff.names,")")
  # res<-do.call(paste,c( as.list(res), sep=""))
  # res
# }
##################################################################
#' @title Compute estimated effects
#' @param Y A numeric named vector where the names are individuals names and the values their phenotype. The names of Y will be matched to the row names of X.
#' @param selec_XXclass A n by mk data.frame of factors
#' with rownames()=individual names, and colnames()=mk selected SNP names
#' additive+dominance: three levels factor
#' female+male+interaction: four levels factor
#'
#' Use function \code{\link{fromeBICtoEstimation}} to get this armument.
#' @param KK a list of one, two or three matrices depending on the models
#'
#' - additive: a n by n matrix, where n=number of individuals, with rownames()=colnames()=individual names
#'
#' - additive+dominance: two n by n matrices, where n=number of individuals, with rownames()=colnames()=individual names
#'
#' - female+male: a n.female by n.female matrix, with rownames()=colnames()=female names and a n.male by n.male matrix, with rownames()=colnames()=male names
#'
#' - female+male+interaction: the same two matrices as the model female+male and a n by n matrix, where n=number of individuals, with rownames()=colnames()=individual names
#' @param cofs A n by q matrix, where n=number of individuals, q=number of fixed effect, with rownames()=individual names and with column names, forbidden head of column names for this matrix "eff1_" and usage of special characters as "*","/","&"
#' @param female A factor of levels female names and length n, only for the last two models
#' @param male A factor of levels male names and length n, only for the last two models
#' @description Estimate the effect of selected SNPs.
#' @return A dataframe with 3 colum: BLUE, Tukey.Class and Frequency. The firt line name is "mu", the names of the other lines are in the form markername_allele.
#' @template examples.donotrun
#' @export
Estimation_allmodels<-function(Y,selec_XXclass,KK,cofs=NULL,female=NULL,male=NULL) {
    #EDIT: adding special character support in marker names ( Olivier Guillaume 2018/07)
    stopifnot(anyDuplicated(colnames(selec_XXclass)) == 0)
    newnames =  gsub("[^a-zA-Z0-9]","_",colnames(selec_XXclass))
    stopifnot(anyDuplicated(newnames) == 0)#different marker names become the same after special character being remplaced
    XX_mrk_names = structure( colnames(selec_XXclass), names = newnames)
    colnames(selec_XXclass) = newnames
    #EDIT END


    stopifnot(!is.null(selec_XXclass))
    nom.effet<-c("eff1","eff2","eff3")
    n<-length(Y)
    nb.effet<-length(KK)
    nb.level.byeffet<-rep(n,nb.effet)  #default for additive and dominance model
    ind<-as.factor(names(Y))  #default for additive and dominance model
    effet<-list()
    for(ki in 1:nb.effet) {
    	effet[[ki]]<-ind
    	}
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
    rm(KK)
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
    #build design matrix
    if(is.null(cofs)){
      design = data.frame(X0,selec_XXclass)
    }else{
      design <- data.frame(X0,as.matrix(cofs),selec_XXclass)
    }
    # model lmekin
    #mod.random<-my.random.lmekin(names(effet))
    	mod.random<-"+(1|ind)"
    xm<-my.model.lmekin(colnames(design),mod.random)


      res.lmekin<-coxme::lmekin(xm, data=data.frame(Y=as.vector(Y), design, effet, ind),
                       varlist = coxme::coxmeMlist(KK.proj) ,method="ML",
                       model = TRUE, x = TRUE, y = TRUE,
                       na.action = na.omit)

      blue = res.lmekin$coefficients$fixed
    #print(xm)
    id<-which(is.na(blue))
    if(length(id)!=0) {
    #this part is due to non estimated levels in the fixed part of the model, this programmation is really not optimised
    	XX<-res.lmekin$x[,-id]
    	res.lmekin<-coxme::lmekin(Y~0+XX+(1|ind), data=data.frame(Y=as.vector(Y), design, effet, ind),
                       varlist = coxme::coxmeMlist(KK.proj) ,method="ML",
                       model = TRUE, x = TRUE, y = TRUE,
                       na.action = na.omit)
    	blue<-res.lmekin$coefficients$fixed
    	names(blue)<-unlist(sapply(names(blue),function(nom){strsplit(nom,"XX")[[1]][2]}))
    	}
    #stopifnot(length(id)==0)
    #change names of blue ... put an  "_" between the name of SNP and the factor level
      #EDIT: the cofactors names were not properly handle (Olivier Guillaume 2018/07)
      #names(blue) = paste0(substr(names(blue),1,nchar(names(blue))-2),"_",substr(names(blue),nchar(names(blue))-1, nchar(names(blue))))
      iNames = ! names(blue) %in% colnames(cofs)
      names(blue)[iNames] = paste0(substr(names(blue)[iNames],1,nchar(names(blue)[iNames])-2),"_",substr(names(blue)[iNames],nchar(names(blue)[iNames])-1, nchar(names(blue)[iNames])))

      names(blue) = gsub("01[|]_10", "_01|10", names(blue))
      names(blue)[1] = "mu"
      varcov.blue = res.lmekin$var
      row.names(varcov.blue) = names(blue)
      colnames(varcov.blue) = names(blue)
      list.snp = colnames(selec_XXclass)
    #
    nddl.res<-length(Y)-length(blue)

    #find  "missing" levels, and compact them
    missing.levels<-unlist(lapply(list.snp,function(s){
    	id<-grep(s,names(blue))
    	s.levels<-levels(as.factor(selec_XXclass[,which(colnames(selec_XXclass)==s)]))
    	names.s.levels<-paste0(s,"_",s.levels)
    	res<-setdiff(names.s.levels,names(blue)[id])
        #compacte levels
    		if(length(res)>1) {
        		lev<-unlist(lapply(res,function(xx){strsplit(xx,paste0(s,"_") )[[1]][2]}))
        		compacte.lev<-do.call(paste,c(lev,list(sep="+")))
        		res<-paste0(s,"_",compacte.lev)
    		}
    	res
    }))

    #complete blue and varblue for "missing" levels, put at 0 as for asreml but with compacted levels
    blue.complete<-c(blue,rep(0,length(missing.levels)))
    names(blue.complete)<-c(names(blue),missing.levels)
    varcov.blue.complete<-cbind(varcov.blue,matrix(0,ncol=length(missing.levels),nrow=nrow(varcov.blue)))
    varcov.blue.complete<-rbind(varcov.blue.complete,matrix(0,nrow=length(missing.levels),ncol=ncol(varcov.blue.complete)))

      #mean comparison of snp genotypes
      res.tukey<-lapply(list.snp,function(i){
        #print(i)
        id<-grep(i,names(blue.complete))
        res = c(NA, names(blue.complete)[id])
        if(length(id)>1){ #do comparison only if more than one level
          mod <- multcomp::parm(blue.complete[id],varcov.blue.complete[id,id],nddl.res)
          contrastes <- multcomp::contrMat(rep(1,length(blue.complete[id]) ),type="Tukey")
          res <- multcomp::glht(model=mod,linfct=contrastes)
          res <- summary(res,test=multcomp::adjusted())
        }
        res
      })
      res.tukey.letter<-unlist(lapply(1:length(res.tukey),function(i){
        if(!is.na(res.tukey[[i]][1])){ # add by Prune 17.06.16
          snp.gen<-names(res.tukey[[i]]$coef)
          mat<-matrix(0,ncol=length(snp.gen),nrow=length(snp.gen) )
          mat[lower.tri(mat)]<-as.vector(res.tukey[[i]]$test$pvalues)
          mat<-mat+t(mat)
          diag(mat)<-1
          colnames(mat)<-snp.gen
          rownames(mat)<-snp.gen
          res<-multcompView::multcompLetters(mat)$Letters
        }else{
          res = rep(NA, length(res.tukey[[i]])-1)
          names(res) = res.tukey[[i]][2:length(res.tukey[[i]])]
        }
        res
      }))
    #change select_XXclass when levels are compacted
    for(s in list.snp){
    	ic<-which(names(selec_XXclass)==s)
    	nb.levels.obs<-length(levels(selec_XXclass[,ic]))
    	nb.levels.compact<-length(grep(s,names(blue.complete)))
    		if(nb.levels.obs!=nb.levels.compact){
    		levels.obs<-levels(selec_XXclass[,ic])
    		levels.compact<-unlist(lapply(names(blue.complete)[grep(s,names(blue.complete))],function(xx){
    			strsplit(xx,paste0(s,"_"))[[1]][2]
    			}))
    		lev.tochange<-setdiff(levels(selec_XXclass[,ic]),levels.compact)
    		lev.compact<-setdiff(levels.compact,levels(selec_XXclass[,ic]))
    		selec_XXclass[,ic]<-as.vector(selec_XXclass[,ic])
    			for(il in 1:length(lev.tochange)){
    				ir<-which(selec_XXclass[,ic]==lev.tochange[il])
    				selec_XXclass[ir,ic]<-lev.compact
    			}
    		}
    }
      #frequencies of classes
      freq<-lapply(selec_XXclass,function(x){
        table(x)/sum(table(x))})
        names.freq<-unlist(lapply(1:length(freq),function(ii){ #same name as blue.complete
        paste0( names(freq)[ii],"_",names(freq[[ii]]) )
      }))
      freq<-unlist(freq)
      names(freq)<-names.freq
      ##########
      RES<-merge(blue.complete,res.tukey.letter,by.x=0,by.y=0,all.x=TRUE)
      names<-RES[,1]
      RES<-RES[,-1]
      rownames(RES)<-names
      RES<-merge(RES,freq,by.x=0,by.y=0,all.x=TRUE)
      names<-RES[,1]
      RES<-RES[,-1]
      rownames(RES)<-names
      colnames(RES)<-c("BLUE","Tukey.Class","Frequency")


    #EDIT: adding special character support in marker names ( Olivier Guillaume 2018/07)
    i = ! rownames(RES) %in% c("mu", colnames(cofs))
    rownames(RES)[i] = paste0(XX_mrk_names[sub("^(.+)_(.+)$","\\1",rownames(RES)[i])],sub("^(.+)_(.+)$","_\\2",rownames(RES)[i]))
    #EDIT END

    #EDIT: reordering the markers so the order will be the same as in selec_XXclass (Olivier Guillaume 2018/07)
    RES = RES[order(match(sub("^(.+)_(.+)$","\\1",rownames(RES)), colnames(selec_XXclass)), na.last = FALSE),]

    RES
}
