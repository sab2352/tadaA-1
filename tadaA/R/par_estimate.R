#' @title Parameter Estimation.
#' @description Joint estimation.
#' @param data the [base_info] returned from [TADA_A_reading_in_annotations], which contains all the allele specific data across all studies
#' @param selected_annotations a vector indicating non-coding annotations whose RRs need to be estimated. e.g., c(2,3) means that the 2nd and 3rd annotations in the [noncoding_annotations] argument of [TADA_A_reading_in_annotations] will have their RRs estimated.
#' @param gene_prior_file a file that has the prior probability of each gene as a risk gene. e.g., "../data/Example_gene_prior.txt".
#' @param optimization_iteration the number of iterations that optim() will perform to estimate RRs.
#' @init initial values for parameter estimation.
#' @hessian hessian in optim(). True or False.
#' @return NULL
#' @export
#' @examples NULL

par_estimate = function(data, selected_annotations, gene_prior_file,optimization_iteration,init,hessian){


	gene_prior = fread(gene_prior_file)


	logP_Zg0 = sumall0(data)

	fr_pi <- function(par){

		    all_rr = par[1:length(selected_annotations)]
	      ###gene_prior$prior <- rep(par[length(par)],nrow(gene_prior))

	      logP_Zg1 = sumall1(data,selected_annotations,all_rr)

	      logP_table<-data.table(logP_Zg1 = logP_Zg1, logP_Zg0 = logP_Zg0, genename = names(data))
		    logP_table <- logP_table[gene_prior, on = "genename"]
		    logP_table <- logP_table[is.finite(logP_table$logP_Zg1) & is.finite(logP_table$logP_Zg0),]

		    idx = match(unique(logP_table$genename),logP_table$genename)
		    idx2 = c(idx,nrow(logP_table) + 1)
			  pr = logP_table[idx,]$prior
			  ll_sum1 <- ll_sum(idx2,logP_table$logP_Zg1,logP_table$logP_Zg0,pr)
			  ll_sum1
			        }

	tmm = proc.time()
	mle0 = optim(init, fr_pi,method="BFGS" ,control=list("fnscale"=-1, "maxit" = optimization_iteration), hessian = hessian)
	tmm2 = proc.time() - tmm
	print(paste0("optim: ",tmm2[3],"s."))
	par_estimate <- c(mle0$par)
	if(hessian){
	fisher_info <- solve(-mle0$hessian)
	prop_sigma <- sqrt(c(diag(fisher_info)))
	upper<-par_estimate+1.96*prop_sigma
  lower<-par_estimate-1.96*prop_sigma
	data.frame(par_estimate,upper,lower)
	    }else{
        data.frame(par_estimate)
	}
}
