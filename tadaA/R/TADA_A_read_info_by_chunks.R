#' @title Read Information through multiple chunks
#' @description A function that is used in the allele-specific model. Read in DNM data, annotation data and adjusted mutation rate data.
#' @param mut_files a vector of files with DNM infomation in a txt format. The first three columns are chromosome, 0-based start and 1-based end, followed by two columns of ref and alt alleles.
#' @param window_file The first 4 columns must be chr, start, end, and site_index of genomic windows. The rest of the columns are features that might affect baseline background mutation rates and that need to be adjusted for.
#' @param mutrate_scaling_files a vector of files that have the scaling factor for each genomic interval in [window_file]. 1st column is site_index, 2nd column is scaling factor of mutation rate. each mutrate_scaling_file matches to one pair of window_file and mut_file. The order of files in [mutrate_scaling_file] should match that in [mut_files].
#' @param sample_sizes a vector of the number of individuals in each study. The order of the numbers should match the order of mutation files in [mut_files].
#' @param gene_prior_file a file that has prior (derived from posterior and prior)for a gene as a risk gene.
#' @param nonAS_noncoding_annotations a vector of non-allele-specific non-coding annotations. Each element is a name of the file that has one non-coding annotation in BED format. Non-coding annotations that are not overlapped with regions in [window_file] will not be used in model fitting.
#' @param AS_noncoding_annotations NA or a list of vectors of allele specific annotations. i.e., Each type of noncoding annotation is an element in the list. An element is comprised of 4 different bed files, corresponding to noncoding annotatins of this type based on the alternative allele, A, T, C, or G. e.g., "spidex_lower10pct_alt_A.bed" is a bed file that has genomic intervals representing the union of all bases which, if mutated to an A allele, have a spidex score lower than 10pct of all possible spidex scores.
#' @param report_proportion Choose the top X percent TADA genes to estimate RR.
#' @param chunk Default is 1. Split into some different chunks, then for each chunk split by genes, then by feature configuration. If split by gene at the first level, implementation would take too long.
#' @param node_n the number of nodes used to run a certain chunk of the code, default is 6.
#' @param mutrate_ref_files a vector of mutrate files in the bigwiggle format. These files have base-level mutation rates to a specific allele, A, C, G, T.
#' @param MPI the index that will add to temp files, useful when running multipe processes at one time. Prefix for temporary files that will be deleted at the end of the pipeline.
#' @return NULL
#' @export
#' @examples NULL
TADA_A_read_info_by_chunks = function(mut_files,
     window_file,
     mutrate_scaling_files,
     sample_sizes,
     gene_prior_file,
     nonAS_noncoding_annotations,
     AS_noncoding_annotations,
     report_proportion,
     chunk ,
     node_n,
     mutrate_ref_files,
     MPI = 1){

    options(warn=-1)
    if(node_n != 1){
    cl <- makeCluster(node_n)
    registerDoParallel(cl)

    f = fread(window_file)
    f = f[order(f$genename),]

    gene = unique(f$genename)

    gn = length(gene)

    step = floor(gn / chunk)

    res = gn - chunk * step

    lst = seq(0,gn,step )

    if(res > 0 ) {
        lst[length(lst)] = gn
        }
    prefix <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
    data_partition = foreach(i = 1:(length(lst)-1), .combine = c ,.packages=c("tadaA")) %dopar% {

       sg = gene[lst[i] + 1]
       eg = gene[lst[i+1]]
       interval_s = min(which(f$genename == sg))
       interval_e = max(which(f$genename == eg))
       fg = f[interval_s:interval_e,]

       fwrite(fg,paste0("tmp/",prefix,"_",i,"_window.bed"),sep="\t")

       TADA_A_read_info(mut_files = mut_files,
                       window_file = paste0("tmp/",prefix,"_", i,"_window.bed"),
                       mutrate_scaling_files =  mutrate_scaling_files,
                       sample_sizes = sample_sizes,
                       gene_prior_file = gene_prior_file,
                       nonAS_noncoding_annotations = nonAS_noncoding_annotations,
                       AS_noncoding_annotations = AS_noncoding_annotations,
                       report_proportion = report_proportion,
                       chunk_partition_num = 1,
                       node_n = node_n,
                       mutrate_ref_files = mutrate_ref_files,
                       MPI = MPI)
        }

    stopCluster(cl)

    }else{

    data_partition = TADA_A_read_info(mut_files = mut_files,
                       window_file = window_file,
                       mutrate_scaling_files =  mutrate_scaling_files,
                       sample_sizes = sample_sizes,
                       gene_prior_file = gene_prior_file,
                       nonAS_noncoding_annotations = nonAS_noncoding_annotations,
                       AS_noncoding_annotations = AS_noncoding_annotations,
                       report_proportion = report_proportion,
                       chunk_partition_num = chunk,
                       node_n = node_n,
                       mutrate_ref_files = mutrate_ref_files,
                       MPI = MPI)
   }

#      saveRDS(tmp,paste0("tmp/",i,".rds"))
#      rm(tmp)
#      return(paste0("chunk",i," finished!"))
#         }

   return(list(base_info = data_partition))
}
