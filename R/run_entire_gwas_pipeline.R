#' @title MLMM, model selection and effects estimation
#' @description
#' Internaly run functions of the mlmm.gwas package:
#' \itemize{
#'   \item \code{\link{mlmm_allmodels}} (GWAS)
#'   \item \code{\link{frommlmm_toebic}}
#'   \item \code{\link{eBIC_allmodels}} (model selection)
#'   \item \code{\link{fromeBICtoEstimation}}
#'   \item \code{\link{Estimation_allmodels}} (effects estimation)
#' }
#' @inheritParams mlmm_allmodels
#' @return
#' A named list with 2 or 3 elements:
#' \itemize{
#'     \item pval: the return value of \code{\link{mlmm_allmodels}}
#'     \item eBic: the return value of \code{\link{eBIC_allmodels}}
#'     \item effects: the return value of \code{\link{Estimation_allmodels}}, only if there is at least one marker in the model selected by lowest eBIC.
#' }
#' @examples
#' data("mlmm.gwas.AD")
#' results <- run_entire_gwas_pipeline(floweringDateAD, list(Xa), list(K.add))
#' @export
run_entire_gwas_pipeline = function(Y, XX, KK, nbchunks=2, maxsteps=20, cofs=NULL, female=NULL, male=NULL, threshold=NULL){
    res = list()

    res$pval = mlmm_allmodels(Y, XX, KK, nbchunks, maxsteps, cofs, female, male)
    sel.XX = frommlmm_toebic(XX, res$pval)
    nb.snp = ncol(XX[[1]])
    res$eBic = eBIC_allmodels(Y, sel.XX, KK, nb.snp, cofs, female, male)
    res$threshold = threshold_allmodels(threshold=threshold, res_mlmm=res$pval)

    if(which.min(res$eBic[, grep("eBIC_", colnames(res$eBic))[1]]) == 1 & length(res$threshold)==1){#nb. the name of the column of the eBIC is variable, ie. eBIC_[lambda]
        warning("The model with the lowest eBic is the null model (without markers used as fixed effects).")
    }else{
        sel.XXclass = fromeBICtoEstimation(sel.XX, res$eBic,res$threshold)
        res$effects = Estimation_allmodels(Y, sel.XXclass, KK, cofs, female, male)
    }
    res
}


# selec_XX = sel.XX
# nb.snp=500