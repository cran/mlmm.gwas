###################################################################
#' @title Selection of the significant marker at each mlmm step according to a threshold
#' @param threshold a value to declare the significant p value. The default value is Bonferroni 0.05
#' @param res.mlmm a list of p-value for each mlmm step  
#' @return A matrix with a line for each mlmm step with p-value of a SNP under the threshold : SNP, p-value, step

threshold_allmodels <- function(threshold=NULL, res_mlmm){

if(is.null(threshold)){threshold=0.05/length(res_mlmm[[2]])}# default threshold

if(length(which(res_mlmm[[2]] < threshold))==0){print("No p-value below the threshold for this trait")}else{

res=matrix(NA,0,3)# build the results matrix
colnames(res)=c("SNP","p-value","MLMM_Step")
res

for(i in 2:length(res_mlmm)){# identify the smallest p-value at each step
# i=2
if(i==2){
if((min(res_mlmm[[i]])<threshold)==FALSE){print(paste0("No p-value below the threshold at step ",(i-1)));break}else{
tab = res_mlmm[[i]][which(res_mlmm[[i]] == min(res_mlmm[[i]]))]
res=rbind(res,c(names(tab),as.vector(tab),(i-1)))}
}else{
if((min(res_mlmm[[i]][-c(1:(i-1))])<threshold)==FALSE){print(paste0("No p-value below the threshold at step ",(i-1)));break}else{
tab = res_mlmm[[i]][-c(1:(i-1))][which(res_mlmm[[i]][-c(1:(i-1))] == min(res_mlmm[[i]][-c(1:(i-1))]))]
res=rbind(res,c(names(tab),as.vector(tab),(i-1)))}
}}

res=as.data.frame(res)
res
}
}

