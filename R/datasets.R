#' Dataset for examples with the mlmm.gwas package, additive and additive+dominance models
#'
#' \itemize{
#'    \item Species: Helianthus annuus
#'    \item Individuals: 125
#'    \item Markers: 500
#' }
#'
#' Variables:
#' \itemize{
#'   \item floweringDateAD: flowering dates in °C.day, time since sowing.
#'   \item Xa: genotype matrix (additive)
#'   \item Xd: genotype matrix (dominance)
#'   \item K.add: "kinship" matrix (additive)
#'   \item K.dom: "kinship" matrix (dominance)
#' }
#'
#'
#' @format
#'
#' \itemize{
#'   \item floweringDateAD: a named numeric of length 444.
#'   \item Xa: a 110x500 numeric matrix
#'   \item Xd: a 110x500 numeric matrix
#'   \item K.add: a 110x110 numeric matrix
#'   \item K.dom: a 110x110 numeric matrix
#' }
#'
#' @docType data
#' @keywords dataset
#' @name mlmm.gwas.AD
#' @usage data(mlmm.gwas.AD)
#' @aliases floweringDateAD Xa Xd K.add K.dom
#' @source \href{https://link.springer.com/article/10.1007/s00122-017-3003-4}{Bonnafous & al. (2017)}
NULL


#' Dataset for examples with the mlmm.gwas package, male+female and male+female+interaction models
#'
#' \itemize{
#'    \item Species: Helianthus annuus
#'    \item Individuals: 303
#'    \item Markers: 500
#' }
#'
#' Variables:
#' \itemize{
#'   \item floweringDateFMI: flowering dates in °C.day, time since sowing.
#'   \item female: names of the female parent of the individuals
#'   \item male: names of the male parent of the individuals
#'   \item hybrid: names of the hybrids (name of female and male)
#'   \item Xf: female genotype matrix (additive)
#'   \item Xm: male genotype matrix (dominance)
#'   \item Xfm: female-male interaction genotype matrix (dominance)
#'   \item K.female: female "kinship" matrix (additive)
#'   \item K.male: male "kinship" matrix (dominance)
#'   \item K.hybrid: hybrid "kinship" matrix (dominance)
#' }
#'
#'
#' @format
#'
#' \itemize{
#'   \item floweringDateFMI: a named numeric vector of length 303.
#'   \item female: a factor of length 303
#'   \item male: a factor of length 303
#'   \item hybrid: a factor of length 303
#'   \item Xf: a 303x500 numeric matrix
#'   \item Xm: a 303x500 numeric matrix
#'   \item Xfm: a 303x500 numeric matrix
#'   \item K.female: 36x36 numeric matrix
#'   \item K.male: 36x36 numeric matrix
#'   \item K.hybrid: 303x303 numeric matrix
#' }
#'
#' @docType data
#' @keywords dataset
#' @name mlmm.gwas.FMI
#' @usage data(mlmm.gwas.FMI)
#' @aliases floweringDateFMI female male hybrid Xf Xm Xfm K.female K.male K.hybrid
#' @source \href{https://link.springer.com/article/10.1007/s00122-017-3003-4}{Bonnafous & al. (2017)}
NULL