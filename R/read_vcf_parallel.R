#' Read VCF: parallel
#' 
#' Read a VCF file across 1 or more threads in parallel.
#' @inheritParams read_vcf
#' @inheritParams GenomicRanges::tileGenome
#' @inheritParams BiocParallel::register
#' @inheritParams downloader
#' @inheritParams download_vcf
#' @inheritParams import_sumstats
#' 
#' @keywords internal
#' @importFrom VariantAnnotation ScanVcfParam VcfFile readVcf rbind
#' @importFrom GenomicRanges seqinfo tileGenome 
#' @importFrom GenomeInfoDb keepSeqlevels 
#' @importFrom data.table rbindlist
#' @returns VCF file.
#' @source 
#' \code{
#' path <- "https://gwas.mrcieu.ac.uk/files/ieu-a-298/ieu-a-298.vcf.gz"
#' 
#' #### Single-threaded ####
#' vcf <- MungeSumstats:::read_vcf_parallel(path = path)
#' #### Parallel ####
#' vcf2 <- MungeSumstats:::read_vcf_parallel(path = path, nThread=11)
#' }
read_vcf_parallel <- function(path,
                              samples = 1,
                              which = NULL,
                              use_params = TRUE,
                              as_datatable = TRUE, 
                              sampled_rows = 1e4L,
                              tilewidth = 1e7L,
                              include_xy = FALSE,
                              
                              download = TRUE,
                              vcf_dir = tempdir(),
                              download_method = "download.file",
                              force_new = FALSE,
                              
                              mt_thresh = 1e8L,
                              nThread = 1,
                              ntile = nThread,
                              verbose = TRUE){
    # echoverseTemplate:::source_all()
    # echoverseTemplate:::args2vars(read_vcf_parallel) 
    requireNamespace("GenomicFiles")
    requireNamespace("VariantAnnotation") 
    
    data.table::setDTthreads(threads = nThread)
    # data.table::getDTthreads()
    ##### Download the VCF to temp ####
    ## This makes importing much faster, 
    ## especially in parallel since you would be makings hundreds 
    ## network calls to OpenGWAS.
    if(nThread>1 && isFALSE(download)){
        messager("Warning: download arg must be TRUE when nThread>1",
                 "to avoid making too many queries to remote file.",
                 "Setting download=TRUE.",
                 v=verbose)
        download <- TRUE
    }
    if(isTRUE(download)){
        vcf_paths <- download_vcf(vcf_url = path,
                                  vcf_dir = vcf_dir,
                                  download_method = download_method,
                                  force_new = force_new,
                                  quiet = !verbose, 
                                  nThread = nThread)
        path <- vcf_paths$save_path
    } 
    #### Optimise query #### 
    if((!is.null(which)) || isTRUE(use_params)){
        #### Make sure file is compressed and indexed ####
        ## File must be indexed in order to use param 
        ## (even if only specifying columns) 
        path <- index_vcf(path = path,
                          verbose = verbose)
        param <- select_vcf_fields(path = path, 
                                   which = which, 
                                   samples = samples,
                                   sampled_rows = sampled_rows,
                                   nThread = nThread)
    } else {
        param <- VariantAnnotation::ScanVcfParam()
    }  
    #### Create VcfFile object #### 
    vcf_file <- VariantAnnotation::VcfFile(file = path, 
                                           index = paste0(path,".tbi"))
    ### Read header ####
    header <- VariantAnnotation::scanVcfHeader(file = vcf_file)
    samples <- VariantAnnotation::samples(header)
    ## Get genome build
    genome <- header@header$contig$assembly[[1]] 
    ## Report total variants
    n_variants <- header@header$SAMPLE$TotalVariants
    if(!is.null(n_variants) & isTRUE(verbose)){
        messager("VCF contains:",
                 formatC(as.integer(n_variants),big.mark = ","),"variant(s)",
                 "x",
                 formatC(nrow(header@header$SAMPLE),big.mark = ","),
                 "sample(s)"
                 )
    }  
    #### Make sure multi-threading makes sense given VCF size ####
    if((as.integer(n_variants)<mt_thresh) && (nThread>1)){
        messager("Processing will be more efficient in single-threaded mode",
                 paste0("when nrows<",mt_thresh,"."),
                 "Temporarily setting nThread=1.",
                 v=verbose)
        nThread <- 1
    }
    #### Single-threaded ####
    t1 <- Sys.time()
    if(nThread==1){
        messager("Reading VCF file: single-threaded",v=verbose) 
        vcf <- VariantAnnotation::readVcf(file = path, 
                                          param = param) 
        if(as_datatable){
            vcf <- vcf2df(vcf = vcf, 
                          add_sample_names = length(samples)!=1,
                          verbose = verbose)
        } 
        
    #### Parallel ####
    ## Uses GenomicFiles to parallelise and speed up reading.
    } else {
        messager("Reading VCF file: multi-threaded",
                 paste0("(",nThread," threads)"),v=verbose)  
        register_cores(workers = nThread,
                       progressbar = verbose)
        ## Check which chromosome are available. 
        xy <- if(isTRUE(include_xy)) c("X","Y") else NULL
        possible_chr <- c(
            c(as.character(seq_len(22)),xy),
            paste0("chr",c(as.character(seq_len(22)),xy))
        )
        used_chr <- possible_chr[possible_chr %in% header@reference]
        ## Tile ranges across the genome 
        tiles <-
            GenomicRanges::seqinfo(vcf_file) |>
            GenomeInfoDb::keepSeqlevels(used_chr) |>
            GenomicRanges::tileGenome(cut.last.tile.in.chrom = FALSE,
                                      ntile = ntile)
        # param1 <- param
        # Create mapping function
        MAP <- function(range,
                        file,
                        param,
                        genome="HG19/GRCh37",
                        as_datatable=TRUE,
                        ...) {
            # param <- get(x = "param1", envir = parent.frame(3L))
            param2 <- VariantAnnotation::ScanVcfParam(which = range,
                                                      fixed = param@fixed,
                                                      info = param@info,
                                                      geno = param@geno,
                                                      samples = param@samples)
            vcf <- VariantAnnotation::readVcf(file = file,
                                              genome = genome,
                                              param = param2)
            if(as_datatable){
                return(vcf2df(vcf = vcf,
                              add_sample_names = FALSE,
                              verbose = FALSE))
            } else {return(vcf)}
        }
        ## Parallelised query
        #### reduceByRange ####
        REDUCE <- if(isTRUE(as_datatable)){
            function(x){data.table::rbindlist(l = x, fill=TRUE)}
            
        } else {
            function(x){do.call(VariantAnnotation::rbind, x)}
        }
        system.time({
            vcf <- GenomicFiles::reduceByRange(ranges = tiles,
                                               files = vcf_file$path,
                                               MAP = MAP,
                                               REDUCE = REDUCE,
                                               iterate = TRUE,
                                               ### Args passed to MAP
                                               param = param,
                                               as_datatable = as_datatable,
                                               genome = genome)
        }) 
        vcf <- REDUCE(vcf)
    }   
    if(verbose) methods::show(round(difftime(Sys.time(),t1),1))
    return(vcf)
}
