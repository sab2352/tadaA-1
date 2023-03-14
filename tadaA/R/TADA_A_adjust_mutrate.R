#' @title Mutation rate adjustment
#' @description Function to calibrate background mutation rate for each DNM study.
#' @param mut_file the file with DNM infomation in a BED format. 0-based start and 1-based end.
#' @param window_file the file with genomic windows. Each line represents one window. The first 4 columns must be chr, start, end, and site_index of genomic windows. The rest of the columns are features that might affect baseline background mutation rates and that need to be adjusted for.
#' @param sample_size the number of individuals.
#' @param scale_features a vector of the names of features that need to be scaled, this is recommended to apply to continuous features, such as GC content, to make easier the interpretation of the effect sizes of features.
#' @param scaling_file_name the name of the file that has the scaling factor for each genomic interval in [window_file]. 1st column is site_index, 2nd column is scaling factor of mutation rate.
#' @param mutrate_mode "regular" means the mutation rate of each window is the total sum of mutation rates over all bases. Under this setting, only in very rare cases a window will have mutation rate equal to 0. The output scaling file would remove these windows, as these windows can't be used in the likelihood estimation. "special" means the mutation rate of each window is the rate of  a specific type of mutations. For example, when analyzing WES data, we would want to use syn mutations to adjust for mutation rates. The mutation rate here should accordingly be the rate of synnonymous mutations. Under this situation, more windows might have mutation rate as 0, but we still need to get the scaling factor for these windows, as the rate of other mutation types may not be 0, thus the bases in these windows are informative. If [genes] set to be not "all", then must set [mutrate_mode] to "special". is included in the list of genes specified by [genes] will be used.
#' @param genes could be "all", or a list of genes. If set to be "all", all windows will be used to adjust mutation rates. Otherwise, only windows with a gene name that is included in the list of genes specified by [genes] will be used.
#' @return NULL
#' @export
#' @examples NULL

TADA_A_adjust_mutrate <- function(mut_file, window_file, sample_size, scale_features, scaling_file_name, mutrate_mode = "regular", genes = "all"){
 
    base_ID = genename = mut_count = mutrate = risk = NULL
    # prefix for temporary files that will be deleted at the end of the pipeline
    prefix <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) # prefix for temporary files that will be deleted at the end of the pipeline
    prefix <- paste("tmp/", prefix, sep = "")

    # make a tmp folder for tmp files
    system("mkdir -p tmp")

    command <- paste("sed -n '1!p' ", window_file, "  | awk '{OFS=\"\t\";print $1,$2,$3,$4}' > ",
                   paste(prefix, "_temp_windows.bed", sep = ""), sep = "")

    system(command)
    print("Made a tmp folder for tmp files.")

    tm = proc.time()
    coverage= as.data.table(do_bedtools_coverage(paste0(prefix, '_temp_windows.bed'),mut_file))[,c(1:3,6,7)]

    tm = proc.time() - tm
    print(paste0("Counted mutation over windows. Time consumed: ",round(tm[3],4),"s"))
    coverage$start = coverage$start - 1
    colnames(coverage) <- c("chr","start","end","site_index","mut_count")
    coverage$site_index = as.numeric(coverage$site_index)

    # read in the window file that has feature annotations which might affect mutation rates
    windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    feature_number <- length(scale_features)

    # merge [window] and [coverage] to link mutation count with feature annotations
    coverage <- coverage[windows[,-1:-3], on="site_index"]


    if(genes == "all"){
        target_genes <- unique(coverage$genename)
      }else{
        target_genes <- fread(genes, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        target_genes <- target_genes[[1]]
      }

 if(mutrate_mode == "regular"  ){
    print("Chose mutrate mode: regular.")
    coverage <- coverage[mutrate !=0]

    if(length(scale_features) != 0){
      for(i in 1:length(scale_features)){
        coverage[[scale_features[i]]] <- as.vector(scale(coverage[[scale_features[i]]]))
      }
    }

    # write the formula
    f <- paste("mut_count ~ ", paste(colnames(coverage)[c(10:(10 + feature_number -1))], collapse = " + "), " + offset(log(2*mutrate*sample_size))", sep = "")
    #fit mutation rate model using a fixed set of features not including histone modification marks.
    tm = proc.time()
    out.offset <- glm(as.formula(f), family = poisson(link="log"), data = coverage)
    scaling <- data.table(site_index = coverage$site_index, scaling_factor = out.offset$fitted.values / (2 * coverage$mutrate * sample_size))
    tm = proc.time() - tm
    print(paste0("Finished fiting mutation rate model and calculating scaling factors. Time consumed: ",round(tm[3],4),"s"))
     } else if(mutrate_mode == "special"){
     print("Chose mutrate mode: special.")
     
     if(!is.na(scale_features)){
      for(i in 1:length(scale_features)){
        coverage[[scale_features[i]]] <- as.vector(scale(coverage[[scale_features[i]]]))
      }
    }
    coverage_mutrate_gt0 <- coverage[mutrate !=0]
    # only keep windows that are assigned to [target_genes]
    coverage_mutrate_gt0 <- coverage_mutrate_gt0[is.element(genename, target_genes)]
    f <- paste("mut_count ~ ", paste(colnames(coverage_mutrate_gt0)[10: (10+ feature_number -1)], collapse = " + "), " + offset(log(2*mutrate*sample_size))", sep = "")
    tm = proc.time()
    out.offset <- glm(as.formula(f), family = poisson(link="log"), data = coverage_mutrate_gt0)
    scaling_factor <- exp(as.matrix(cbind(1, coverage[,10 : (10 + feature_number -1)])) %*% out.offset$coefficients)
    scaling <- data.table(site_index = coverage$site_index, scaling_factor = scaling_factor)
    tm = proc.time() - tm
    print(paste0("Finished fiting mutation rate model and calculating scaling factors. Time consumed: ",round(tm[3],4),"s"))
    colnames(scaling) <- c("site_index", "scaling_factor")
  }


  fwrite(scaling, scaling_file_name, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    # remove intermediate files
  system(paste("rm ", prefix, "_temp*", sep = ""))
  print(paste0("Removed temporary files."))
  # the return value is the effect sizes of feature annotations on background mutation rates.
  return(summary(out.offset)$coeff)
}
