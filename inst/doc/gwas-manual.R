## ----setup, include = FALSE----------------------------------------------
    knitr::opts_chunk$set(
      collapse = TRUE,
      comment = "#>"
    )

## ----data----------------------------------------------------------------
    library(mlmm.gwas)
    data("mlmm.gwas.AD")
    ls()

## ----matrixFormat, include=FALSE-----------------------------------------
    # Change le formatage par d??faut des matrices
    # On le met ici pour pas qu'il soit visible quand on fait ls() ci dessus
    knit_print.matrix <- function(x, ...){
        knitr::knit_print(as.data.frame(x), ...)
    }

## ----floweringDateAD-----------------------------------------------------
    head(floweringDateAD)

## ----X-------------------------------------------------------------------
    Xa[1:5,1:5]

## ----K-------------------------------------------------------------------
    ##You can get K from X with the following command line
    ##We are not doing it here because the X matrix included in data("mlmm.gwas.AD")
    ##is a subset of a larger matrix.
    #K.add = Xa %*% t(Xa)

    K.add[1:5,1:5]

## ----mlmm----------------------------------------------------------------
    res_mlmm = mlmm_allmodels(floweringDateAD, list(Xa), list(K.add))

## ----manhattan, fig.height=5, fig.width=5--------------------------------
    manhattan.plot( res_mlmm )

## ----eBIC----------------------------------------------------------------
    sel_XX = frommlmm_toebic(list(Xa), res_mlmm)
    res_eBIC = eBIC_allmodels(floweringDateAD, sel_XX, list(K.add), ncol(Xa))
    
    res_eBIC

## ----threshold-----------------------------------------------------------
    res_threshold <- threshold_allmodels(threshold=NULL, res_mlmm)

## ----estimation----------------------------------------------------------
    sel_XXclass = fromeBICtoEstimation(sel_XX, res_eBIC, res_threshold)
    effects = Estimation_allmodels(floweringDateAD, sel_XXclass, list(K.add))
    
    effects

## ----boxplot-------------------------------------------------------------
    genotypes.boxplot(Xa, floweringDateAD, "SNP303", effects)

## ----completeExample, eval=FALSE-----------------------------------------
#      library(mlmm.gwas)
#      data("mlmm.gwas.AD")
#  
#      res_mlmm = mlmm_allmodels(floweringDateAD, list(Xa), list(K.add))
#  
#      manhattan.plot( res_mlmm )
#  
#      sel_XX = frommlmm_toebic(list(Xa), res_mlmm)
#      res_eBIC = eBIC_allmodels(floweringDateAD, sel_XX, list(K.add), ncol(Xa))
#  
#      res_threshold <- threshold_allmodels(threshold=NULL, res_mlmm)
#  
#      sel_XXclass = fromeBICtoEstimation(sel_XX, res_eBIC, res_threshold)
#      effects = Estimation_allmodels(floweringDateAD, sel_XXclass, list(K.add))
#  
#      genotypes.boxplot(Xa, floweringDateAD, "SNP303", effects)

## ----runEntirePipeline---------------------------------------------------
    results = run_entire_gwas_pipeline(floweringDateAD, list(Xa), list(K.add), threshold=NULL)
    names(results)

