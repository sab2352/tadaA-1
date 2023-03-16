#' @title Read Information
#' @description A function that is used in the allele-specific model. Read in DNM data, annotation data and adjusted mutation rate data.
#' @param mut_files a vector of files with DNM infomation in a txt format. The first three columns are chromosome, 0-based start and 1-based end, followed by two columns of ref and alt alleles. 
#' @param window_file The first 4 columns must be chr, start, end, and site_index of genomic windows. The rest of the columns are features that might affect baseline background mutation rates and that need to be adjusted for.
#' @param mutrate_scaling_files a vector of files that have the scaling factor for each genomic interval in [window_file]. 1st column is site_index, 2nd column is scaling factor of mutation rate. each mutrate_scaling_file matches to one pair of window_file and mut_file. The order of files in [mutrate_scaling_file] should match that in [mut_files].
#' @param sample_sizes a vector of the number of individuals in each study. The order of the numbers should match the order of mutation files in [mut_files].
#' @param gene_prior_file a file that has prior (derived from posterior and prior)for a gene as a risk gene.
#' @param nonAS_noncoding_annotations a vector of non-allele-specific non-coding annotations. Each element is a name of the file that has one non-coding annotation in BED format. Non-coding annotations that are not overlapped with regions in [window_file] will not be used in model fitting.
#' @param AS_noncoding_annotations NA or a list of vectors of allele specific annotations. i.e., Each type of noncoding annotation is an element in the list. An element is comprised of 4 different bed files, corresponding to noncoding annotatins of this type based on the alternative allele, A, T, C, or G. e.g., "spidex_lower10pct_alt_A.bed" is a bed file that has genomic intervals representing the union of all bases which, if mutated to an A allele, have a spidex score lower than 10pct of all possible spidex scores.
#' @param report_proportion Choose the top X percent TADA genes to estimate RR.
#' @param chunk_partition_num Default is 1. Split into some different chunks, then for each chunk split by genes, then by feature configuration. If split by gene at the first level, implementation would take too long.
#' @param node_n the number of nodes used to run a certain chunk of the code, default is 6.
#' @param mutrate_ref_files a vector of mutrate files in the bigwiggle format. These files have base-level mutation rates to a specific allele, A, C, G, T.
#' @param MPI the index that will add to temp files, useful when running multipe processes at one time. Prefix for temporary files that will be deleted at the end of the pipeline.
#' @return NULL
#' @export
#' @examples NULL
TADA_A_read_info <- function(mut_files,
                             window_file,
                             mutrate_scaling_files,
                             sample_sizes,
                             gene_prior_file,
                             nonAS_noncoding_annotations,
                             AS_noncoding_annotations,
                             report_proportion,
                             chunk_partition_num = 1,
                             node_n = 6,
                             mutrate_ref_files,
                             MPI = 1){

  V5 = NULL
  prefix <- prefix <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, MPI, sep = "")


  system("mkdir -p tmp")
  print("Made a tmp folder for tmp files.")

  #mut <- fread(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(windows) = c('chr','start','end','site_index','genename','mutrate','coding','promoter','GC_content','div_score')
  coverage <- windows
  print(paste0('Read in window files.'))

  # the number of genomic windows in [mutrate_scaling] is
  # less than the number of windows in [windows] because there are
  # a few windows with mutration rate equal to 0, and thus removed.
  for(i in 1:length(mut_files)){
    mutrate_scaling <- fread(mutrate_scaling_files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    #system(paste("echo \"Finished reading in mutrate scaling file ", mutrate_scaling_files[i], ".\"", sep = ""))
    #system("date")
    print(paste0("Read in mutrate scaling file ", mutrate_scaling_files[i]))
    coverage <- coverage[mutrate_scaling, on = "site_index"]
    coverage <- coverage[complete.cases(coverage)] # release memory
    colnames(coverage)[length(colnames(coverage))] <- paste("scaling_factor_study_", i, sep = "")
    rm(mutrate_scaling) # release memory
  }

  gene_prior <- fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  print(paste0("Read in gene prior file."))

  # merge gene prior info
  coverage <- coverage[gene_prior, on = "genename"]
  coverage <-coverage[complete.cases(coverage)]

  print(paste0("Chose the top ",round(report_proportion*100,2),"% genes for parameter estimation."))

  # select genes based on TADA prior probability and [report_proportion]
  if(report_proportion !=1){
    genes_for_report <- gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report <- genes_for_report[1:floor(nrow(genes_for_report)*report_proportion)]
    coverage <- coverage[genes_for_report, on = "genename"]
    coverage <-coverage[complete.cases(coverage)]
  }else{
    genes_for_report  <- gene_prior[,1] # choose all the genes in TADA coding table
  }


  # now need to extropolate window mutation file to base level mutation file
  coverage <-coverage[,c("chr","start","end",paste("scaling_factor_study_", seq(1,length(mut_files)), sep = ""),"genename"),with = FALSE]
  coverage$ID <- paste(coverage$genename, coverage$start, sep = "_")

  total_rows <- nrow(coverage)
  interval <- floor(total_rows/chunk_partition_num)
  data_bins <- c(rep(seq(1,chunk_partition_num), each = interval),rep(chunk_partition_num, total_rows -interval*chunk_partition_num))

  # split into 20 different chunks, then for each chunk split by genes, then by feature configuration. If split by gene at the first level, implementation would take too long
  coverage <- split(coverage, data_bins)

  #funtion to expand windows to bases
  window_expansion <- function(table_row){
    start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
    data.frame(table_row[1], start, start+1, paste(table_row["chr"],table_row["genename"],start,sep = "_"), table_row["ID"])
  }

  options(warn=-1)
  if(node_n != 1){
    cl <- makeCluster(node_n)
  }
  # use parallel computing and rbindlist to increase the speed by thousands of times.
  environment(window_expansion) <- .GlobalEnv

  # get nonAS feature number
  if(is.na(nonAS_noncoding_annotations[1])){
    nonAS_feature_number <- 0
  }else{
    nonAS_feature_number <- length(nonAS_noncoding_annotations)
  }

  # get AS feature number
  if(is.na(AS_noncoding_annotations[1])){
    AS_feature_number <- 0
  }else{
    AS_feature_number <- length(AS_noncoding_annotations)
  }

  # get total feature number
  feature_number = nonAS_feature_number + AS_feature_number
  # function to get effective information of each element of partition_by_gene
  # These information are those necessary to compute log-likelihood in the optimization function
  partition_feature <- function(pbg){
    # input is one element of the list of partition_by_gene
    pbg_split <- split(pbg, pbg[,4:(4 + feature_number - 1)],drop = TRUE)
    feature_combination_number <- length(pbg_split)
    # this function below is different from the function used in dealing with dataset without reading by chunk. 
      # Here, prior is not incoporated at this step.
    info_for_each_feature <- function(feature_set){
      list(feature_vector = c(as.numeric(feature_set[1,4:(4 + feature_number - 1)])), 
           sum_mut_rate_count = sum(feature_set$mut_count*log(feature_set$adjusted_base_mutrate)), 
           sum_mut_rate = sum(feature_set$adjusted_base_mutrate), 
           sum_mut_count = sum(feature_set$mut_count), 
           log_fcount = sum(log(factorial(feature_set$mut_count))))
    }
    sapply(pbg_split, info_for_each_feature,simplify = FALSE)
  }

  # build a list to store data
  data_partition <-list()
  alt_letters <- c("A","C","G","T")

  # here only deals with the noncoding parts that are within 10kb of TSSs of genes. Will deal with
  for(i in 1:chunk_partition_num){

    tm = proc.time()
    # split coverage_noncoding to 10 parts, and for each part, generate feature table (which will be used for )
    if(node_n != 1){
      coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage[[i]], 1, window_expansion))
    }else{
      coverage_noncoding_for_base_mutrate <- rbindlist(apply(coverage[[i]], 1, window_expansion))
    }

    #system(paste("echo \"Finished partitioning base-level coordinates data at Round ", i, ".\"", sep = ""))
    #system("date")
    tm = proc.time() - tm
    print(paste0("Patitioned genomic window to base-level coordinates at Round ", i,
                 ". Time consumed: ",round(tm[3],4),"s."))

    colnames(coverage_noncoding_for_base_mutrate) <- c("chr","start","end","base_ID","ID")
    coverage_noncoding_for_base_mutrate$start <- as.integer(coverage_noncoding_for_base_mutrate$start)
    coverage_noncoding_for_base_mutrate$end <- as.integer(coverage_noncoding_for_base_mutrate$end)
	
    # write out a bed file to get base-level mutation rates
    fwrite(coverage_noncoding_for_base_mutrate[,1:4],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    # read in allele-specific base-level mutation rate
      
    tm = proc.time()
    for(j in 1:length(mutrate_ref_files)){
      command <- paste("bigWigAverageOverBed ",
                       mutrate_ref_files[j], " ", paste(prefix, 
                                                        "_temp_for_mutrate.bed", sep = ""), " ", 
                       paste(prefix, "_temp_for_mutrate.bed.mutrate.txtttt", sep = "" ), sep = "")
      system(command)

      #######################  
        
      command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
      system(command)
      coverage_noncoding_base_mutrate <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
      colnames(coverage_noncoding_base_mutrate) <- c("base_ID",paste("base_mutrate_alt_", alt_letters[j], sep = ""))
      coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[coverage_noncoding_base_mutrate, on = "base_ID"]
    }
    tm = proc.time() - tm
    print(paste0("Obtained un-calibrated base-level mutation rate for alt alleles."," Time consumed: ",round(tm[3],4),".s"))

    # read in non allele-specific epigenomic annotations
    epi_ID = 1
    if (!is.na(nonAS_noncoding_annotations)[1]){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
      for(epi in nonAS_noncoding_annotations){
        tm = proc.time()
        base_in_epi = as.data.table(do_bedtools_coverage(paste0(prefix, "_temp_for_mutrate.bed"),epi))[,c(6,7)]
        colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, sep = "_"))
        coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[base_in_epi, on = "base_ID"]
        tm = proc.time() - tm
        print(paste0("Added non-allele specific annotations ", epi_ID ,"."," Time consumed: ",round(tm[3],4),"s."))
        #system(paste("echo \"Read in non-allele specific annotations ", epi_ID, ".\"", sep = ""))
        #system("date")
        epi_ID <- epi_ID + 1
      }
    }

    # read in allele-specific epigenomic annotations, now the base assignment for such annotations is based on the 50-bp windows.
    # Should not be a big problem as distal introns are assigned based on refseq_ID for intron regions, should be consistent with assignment by spidex.
    # though there is a small chance that a gene's distal intron, which is close to the promoter of another gene and within 10 kb from the TSSs of the latter gene, may be mis-assigned to the latter gene.
    # To completely correct for this issue, need to allow epigenetic annotation to take its own gene assignment, 
    # which might be necessary in some situations under strict criteria, such as splicing annotaion.

    if (!is.na(AS_noncoding_annotations)[1]){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
      for(epi in AS_noncoding_annotations){
        tm = proc.time()
        for(k in 1:length(epi)){
          base_in_epi = as.data.table(do_bedtools_coverage(paste0(prefix, "_temp_for_mutrate.bed"),epi[k]))[,c(6,7)]
          colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, alt_letters[k], sep = "_"))
          coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[base_in_epi, on = "base_ID"]
        }
        tm = proc.time() - tm
        print(paste0("Added allele specific annotations ", epi_ID," .Time consumed: ",round(tm[3],4),"s."))
        epi_ID <- epi_ID + 1
      }
    }

    print("Begin collapsing data...")

    
    # now for each study, read in data, and collapse data based on annotation configuration
    for(m in 1:length(mut_files)){
        
      tm = proc.time()
      mut <- fread(mut_files[m], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      print(paste0("Read in mutation file ", mut_files[m], "."))
        
      for(letter in alt_letters){
        if(!is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", 
                                                                                                                    letter, sep = ""), 
                                                                                              paste("Anno_", seq(1, nonAS_feature_number),
                                                                                                    sep = ""), 
                                                                                              paste("Anno_", seq(nonAS_feature_number +1, 
                                                                                                                 nonAS_feature_number + 
                                                                                                                 AS_feature_number), "_", 
                                                                                                    letter ,sep = "")), 
                                                                                              with = FALSE]
            
        }else if(!is.na(nonAS_noncoding_annotations)[1] & is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(1, nonAS_feature_number), sep = "")), with = FALSE]
        }else if(is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }

        # a very important step here is removing bases that have adjusted_mutrate_base 0. This happens when the allele of 
        # the mutrate we are using is just the reference allele.
        # By doing this, we make the computation of likelihood valid, also we automatically removed bases with nonAS annotations 
        # but with mutant allele the same with ref allele at this step
        #tm = proc.time()
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[coverage_noncoding_for_base_mutrate_temp[[paste("base_mutrate_alt_", 
                                                                                                                                             letter, sep = "")]] !=0]
    
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[coverage[[i]][,c(paste("scaling_factor_study_", 
                                                        m, sep = ""), "genename","ID"), with = FALSE],on = "ID"]
        coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate <- coverage_noncoding_for_base_mutrate_temp[, 
                    paste("base_mutrate_alt_", letter, sep = ""), with = FALSE] * coverage_noncoding_for_base_mutrate_temp[, 
                    paste("scaling_factor_study_", m, sep = ""), with = FALSE] * 2 * sample_sizes[m]
        #tm = proc.time() - tm
        #print(paste0("Removed bases that have adjusted_mutrate_base 0 and adjusted base mutrate. Time consumed: ",round(tm[3],4),"s."))
         
        #tm = proc.time()
        if(!is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[, c("base_ID", "adjusted_base_mutrate", 
                                                                                                   "genename", paste("Anno_", 
                                                                                                  seq(1, 
                                                                                                  nonAS_feature_number),
                                                                                                  sep = ""), paste("Anno_", 
                                                                                                  seq(nonAS_feature_number +1, 
                                                                                                  nonAS_feature_number + 
                                                                                                  AS_feature_number), "_", 
                                                                                                                   letter ,sep = "")), 
                                                                                                  with = FALSE]
        }else if(!is.na(nonAS_noncoding_annotations)[1] & is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[, c("base_ID", "adjusted_base_mutrate", "genename", paste("Anno_", seq(1, nonAS_feature_number), sep = "")), with = FALSE]
        }else if(is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[, c("base_ID", "adjusted_base_mutrate", 
                                "genename", paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number),
                                                  "_", letter ,sep = "")), with = FALSE]
        }
        #tm = proc.time() - tm
        #print(paste0("Removed bases with nonAS annotations. Time consumed: ",round(tm[3],4),".s"))

        mut_allele <- mut[V5 == letter]
        # the subsequent two lines of code are used to prevent a outputing bug when using fwrite. (85000000 to be written as 8.5e7)
        mut_allele$V2 <- as.integer(mut_allele$V2)
        mut_allele$V3 <- as.integer(mut_allele$V3)
        fwrite(mut_allele, paste(prefix, "_temp_mut_allele.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
        fwrite(mut_allele[,1:3], paste(prefix, "_temp_mut_allele_ct.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
        ################################ 
        #tm = proc.time()
        base_with_mut = as.data.table(do_bedtools_coverage(paste0(prefix, "_temp_for_mutrate.bed"),
                                                 paste0(prefix, "_temp_mut_allele_ct.bed")))[,6:7]

        colnames(base_with_mut) <- c("base_ID", "mut_count")

        #tm = proc.time() - tm
        #print(paste0("Combined base with mut count. Time consumed: ",round(tm[3],4),".s"))
          
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[base_with_mut, on = "base_ID"]
        
        #coverage_noncoding_for_base_mutrate_temp$Anno_1 = as.numeric(coverage_noncoding_for_base_mutrate_temp$Anno_1)
        #coverage_noncoding_for_base_mutrate_temp$mut_count = as.numeric(coverage_noncoding_for_base_mutrate_temp$mut_count)
        # have to collpase data at this point, otherwise, the I will run out of RAM if process 1000 genes every time.
      
        anno_count <- rep(0, nrow(coverage_noncoding_for_base_mutrate_temp))
        for(p in 1:feature_number){
          anno_count <- anno_count + coverage_noncoding_for_base_mutrate_temp[[(4 + p -1)]]
        }
        # remove rows(bases) that don't have any non-coding features, this could save a lot of RAM, so I could use smaller partition number which would greatly accelerate speed when read in data for all genes.
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[anno_count >0,]

        # first partition by gene for the current chunk
        #tm = proc.time()
        coverage_noncoding_for_base_mutrate_temp<- split(coverage_noncoding_for_base_mutrate_temp, 
                                                         coverage_noncoding_for_base_mutrate_temp$genename)
        #tm = proc.time() - tm
        #print(paste0("Finished partition by gene for the current chunk. Time consumed: ",round(tm[3],4),".s"))
          
        # then partition by feature configuration for each gene in the current chunk
        #tm = proc.time() 
        coverage_noncoding_for_base_mutrate_temp <- sapply(coverage_noncoding_for_base_mutrate_temp, partition_feature, simplify = FALSE)
        #tm = proc.time() - tm
        #print(paste0("Finished partition by feature configuration for each gene in the current chunk. Time consumed: ",round(tm[3],4),".s"))

        # add compact data
        data_partition <- append(data_partition, coverage_noncoding_for_base_mutrate_temp)
        rm(coverage_noncoding_for_base_mutrate_temp) # release memory
        #print(paste0("Finished read in mutation data and make them into the compact format for Study ", m, " and allele ", letter,"."))
      }
    tm = proc.time() - tm
    print(paste0("Collapsed data based on annotations for Study ", m, " .Time consumed: ",round(tm[3],4),".s"))
        
    }
  }

  if(node_n != 1){
    stopCluster(cl)
  }
  options(warn = 0)
  # remove temporary files
 # system(paste("rm ", prefix, "_temp*", sep = ""))
 # print("Temp files cleaned and data recording finished!")
 #system("date")
  #return(list(base_info = data_partition))
  return(data_partition)
}
