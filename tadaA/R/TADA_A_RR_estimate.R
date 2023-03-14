#' TADA-A RR Estimation
#' @description Estimate relative risks.
#' @param data the [base_info] returned from [TADA_A_reading_in_annotations], which contains all the allele specific data across all studies
#' @param selected_annotations a vector indicating non-coding annotations whose RRs need to be estimated. e.g., c(2,3) means that the 2nd and 3rd annotations in the [noncoding_annotations] argument of [TADA_A_reading_in_annotations] will have their RRs estimated.
#' @param gene_prior_file a file that has the prior probability of each gene as a risk gene. e.g., "../data/Example_gene_prior.txt".
#' @param optimization_iteration the number of iterations that optim() will perform to estimate RRs.
#' @param mode "regular", or "single_fast". "single_fast" is used when estimating RR from only one annotation ([data] only recoreded one annotation) of lots of genes (e.g., all genes), would be 5 times faster.
#' @return NULL
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @export
#' @examples NULL
#' @useDynLib(tadaA)


TADA_A_RR_estimate <-function(data, selected_annotations, gene_prior_file, optimization_iteration = 2000, mode = "regular",hessian=TRUE){

# get the piror probability of genes.
tm = proc.time()
gene_prior = fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(gene_prior) = c("genename", "prior")
gene_prior = gene_prior[order(gene_prior$genename),]
gene_prior$prior <- as.numeric(gene_prior$prior)
gene_prior$genename = as.character(gene_prior$genename)
tm = proc.time() - tm
print(paste0("Read in gene prior file. Time consumed: ", round(tm[3],4),"s."))

logP_Zg0 = sumall0(data)

fr <- function(x){
  all_rr = x

  logP_Zg1 = sumall1(data,selected_annotations,all_rr)

  logP_table<-data.table(logP_Zg1 = logP_Zg1, logP_Zg0 = logP_Zg0, genename = names(data))
  logP_table <- logP_table[gene_prior, on = "genename"]
  logP_table <- logP_table[complete.cases(logP_table)]
  idx = match(unique(logP_table$genename),logP_table$genename)
  idx2 = c(idx,nrow(logP_table) + 1)
  pr = logP_table[idx,]$prior
  ll_sum1 <- ll_sum(idx2,logP_table$logP_Zg1,logP_table$logP_Zg0,pr)
  ll_sum1
}

if(length(selected_annotations) == 1){
    tm =  proc.time()
    mle = optim(rep(0.1, 1), fr, method = "Brent", lower = -1, upper = 10,
                                 control=list("fnscale"=-1, "maxit" = optimization_iteration), hessian = hessian)
    par_estimate <- c(mle$par)
if(hessian){
    fisher_info <- solve(-mle$hessian)
    prop_sigma <- sqrt(c(diag(fisher_info)))
    upper<-par_estimate+1.96*prop_sigma
    lower<-par_estimate-1.96*prop_sigma
    res <-data.frame(par_estimate, lower, upper) } else {res <-data.frame(par_estimate) }
    } else {

    res = par_estimate(data=data, selected_annotations=selected_annotations,
                   gene_prior_file = gene_prior_file,
                   optimization_iteration = optimization_iteration,
                   #init = c(rep(.1,length(selected_annotations)), .5),hessian = hessian)
                   init = c(rep(.1,length(selected_annotations))),hessian = hessian)
}

print("Finished parameter estimation!")
return(res)

    }
