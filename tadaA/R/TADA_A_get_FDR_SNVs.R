#' TADA_A_get_FDR_SNVs
#' @description Caculate posterior probabilities after considering SNV annotations
#' @param data the [base_info] returned from [TADA_A_read_info], each element of the list has compact base-level information for each gene.
#' @param selected_annotations determines which features in the [data] are going to be used for calculating bayes factors. e.g., if [selected_annotations] is c(1,3). Then the first and third element of feature_vector in the [base_info] of [data] will be selected.
#' @param gama the estimated relative risks of [selected_annotations]. gama is estimated from [TADA_A_RR_estimate].
#' @gene_prior_file The original prior file
#' @return NULL
#' @import Rcpp
#' @export
#' @examples NULL
#' @useDynLib(tadaA)
#' @importFrom Rcpp sourceCpp

TADA_A_get_FDR_SNVs <- function(data,gama,selected_annotations,gene_prior_file,optimization_iteration=2000){

  gene_prior = fread(gene_prior_file)
  #gene_prior$prior = rep(gama[length(selected_annotations)+1],nrow(gene_prior))
  logP_Zg0 = sumall0(data)
  post_fr <- function(x){
    all_rr = x
    logP_Zg1 = sumall1(data,selected_annotations,all_rr)
    logP_table<-data.table(logP_Zg1 = logP_Zg1, logP_Zg0 = logP_Zg0, genename = names(data))
    logP_table <- logP_table[gene_prior, on = "genename"]
    logP_table <- logP_table[complete.cases(logP_table)]
    u = unique(logP_table$genename)
    idx = match(u,logP_table$genename)
    idx2 = c(idx,nrow(logP_table) + 1)
    pr = logP_table[idx,]$prior
    post = post(idx2,logP_table$logP_Zg1,logP_table$logP_Zg0,pr)
    post.dt = data.table(genename =u,prior = post )
    post.dt
  }

  #g = post_fr(gama[-(length(selected_annotations)+1)])
  g = post_fr(gama)

  g$q0 = 1- g$prior

  g2 = g[order(g$q0),]
  rm(g)

  FDR = c()
  for (i in 1:nrow(g2)) FDR[i] <- sum(g2$q0[1:i]) / i
  g2$FDR = FDR

  list(FDR_ls_0.1=g2[g2$FDR<0.1,],FDR_all=g2)

}
