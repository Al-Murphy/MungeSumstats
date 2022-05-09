#' VCF to DF
#' 
#' Function to convert a \pkg{VariantAnnotation}
#'  \code{CollapsedVCF}/\code{ExpandedVCF}
#'  object to a \code{data.frame}. 
#' @param vcf Variant Call Format (VCF) file imported into R 
#' as a \pkg{VariantAnnotation} 
#' \link[VariantAnnotation]{CollapsedVCF}/
#' \link[VariantAnnotation]{ExpandedVCF} object. 
#' @param add_sample_names Append sample names to column names 
#' (e.g. "EZ" --> "EZ_ubm-a-2929").
#' @param add_rowranges Include \code{rowRanges} from VCF as well.
#' @param drop_empty_cols Drop columns that are filled entirely with: 
#' \code{NA}, \code{"."}, or \code{""}. 
#' @param unique_cols Only keep uniquely named columns.
#' @param unique_rows Only keep unique rows.
#' @param unlist_cols If any columns are lists instead of vectors, unlist them.
#' Required to be \code{TRUE} when \code{unique_rows=TRUE}.
#' @param verbose Print messages.
#' @inheritParams check_empty_cols
#' 
#' @source 
#' \href{https://gist.github.com/zhujack/849b75f5a8305edaeca1001dfb9c3fe9}{
#' Original code source}
#' @source {
#' #### vcfR ####
#' if(!require("pinfsc50")) install.packages("pinfsc50")
#' vcf_file <- system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50")
#' vcf <- read.vcfR( vcf_file, verbose = FALSE )
#' vcf_df_list <- vcfR::vcfR2tidy(vcf, single_frame=TRUE)
#' vcf_df <- data.table::data.table(vcf_df_list$dat)
#' }
#' @return data.frame version of VCF 
#' 
#' @export
#' @importFrom utils type.convert 
#' @importFrom data.table as.data.table
#' @importFrom methods slotNames is show
#' @examples   
#' #### VariantAnnotation ####
#' # path <- "https://github.com/brentp/vcfanno/raw/master/example/exac.vcf.gz"
#' path <- system.file("extdata", "ALSvcf.vcf",
#'                     package = "MungeSumstats")
#'                     
#' vcf <- VariantAnnotation::readVcf(file = path)
#' vcf_df <- MungeSumstats:::vcf2df(vcf = vcf)
vcf2df <- function(vcf, 
                   add_sample_names=TRUE,
                   add_rowranges=TRUE,
                   drop_empty_cols=TRUE,
                   unique_cols=TRUE,
                   unique_rows=TRUE,
                   unlist_cols=TRUE,
                   sampled_rows=NULL,
                   verbose=TRUE) {
    requireNamespace("VariantAnnotation")
    requireNamespace("MatrixGenerics")
  
    messager("Converting VCF to data.table.",v=verbose) 
    #### Expand VCF ####
    if (methods::is(vcf,"CollapsedVCF")) {
        messager('Expanding VCF first, so number of rows may increase.',
                 v=verbose)
        ## Not all VCFs work with this function
        vcf <- tryCatch({
            VariantAnnotation::expand(x = vcf)
        }, error = function(e){message(e); vcf})
    }  
    #### Automatically determine whether to add sample names ####
    if(is.null(add_sample_names)){
        header <- VariantAnnotation::header(x = vcf)
        samples <- VariantAnnotation::samples(header)
        add_sample_names <- length(samples)!=1
    } 
    #### .anncols function ####
    .anncols = function(anncol,headerstring) {
        anncols <- strsplit(sub("Functional annotations: '",'',
                               headerstring),' \\| ')[[1]]
        dfannempty <- data.frame(matrix(vector(), 0, length(anncols),
                                       dimnames=list(c(), anncols)),
                                stringsAsFactors=FALSE)
        yy <- data.frame(
            suppressWarnings(
                do.call(
                    rbind,
                    c(dfannempty,lapply(lapply(anncol,`[`,1),
                                        function(x){strsplit(x,'\\|')[[1]]})
                      )
                    )
                ),
            stringsAsFactors=FALSE)
        yy <- data.frame(lapply(yy,type.convert))
        colnames(yy) <- paste("ANN",anncols,sep="_")
        return(yy)
    }
    #### convert to data.table #### 
    t1 <- Sys.time()
    #### Get rowranges only if available ####
    if("rowRanges" %in% methods::slotNames(vcf)){
        gr <- MatrixGenerics::rowRanges(vcf)
    } else {
        gr <- NULL
    }
    df <- data.table::data.table(
        ID = names(gr),
        if(add_rowranges) granges_to_dt(gr = gr) else NULL,
        DF_to_dt(DF = VariantAnnotation::info(vcf)) 
    )
    remove(gr)
    #### Parse ANN column ####
    if('ANN' %in% colnames(df)) {
        dfann <- .anncols(
            anncol = df$ANN,
            headerstring = VariantAnnotation::info(
                VariantAnnotation::header(vcf)
                )['ANN',]$Description
        )
        df <- df[,colnames(df)!="ANN"]
        df <- cbind(df,dfann)
    }
    ##### Convert geno data ##### 
    ## Works better for VCFHeader
    if(methods::is(vcf,"VCFHeader")){
        tmp <- data.table::as.data.table(VariantAnnotation::geno(vcf), 
                                         keep.rownames = "name")
        # tmp <- DF_to_dt(DF = VariantAnnotation::geno(x))    
    } else {
        ## Works better for VCF
        n <- names(VariantAnnotation::geno(vcf))
        #### Avoid parsing redundant columns ####
        if("ID" %in% colnames(df)) n <- n[n!="ID"]
        tmp <- lapply(n,function(col) {
            coldat <- VariantAnnotation::geno(vcf)[[col]]  
            #### Convert 3D matrix --> 2D ####
            ## This sometimes happens by accident with 
            ## VariantAnnotation::writeVcf
            if(length(dim(coldat))==3){
                rn <- rownames(coldat)
                cn <- rep(colnames(coldat), dim(coldat)[3])
                if(length(unique(colnames(coldat)))==1){
                    coldat <- matrix(data = coldat[,,1, drop=FALSE])
                    rownames(coldat) <- rn
                    colnames(coldat) <- cn[1]
                } else {
                    coldat <- matrix(data = coldat,
                           nrow = dim(coldat)[1], 
                           ncol = dim(coldat)[2] * dim(coldat)[3])
                    rownames(coldat) <- rn
                    colnames(coldat) <- cn
                }  
            }
            #### Drop empty columns within each matrix ####
            if(isTRUE(drop_empty_cols)){
                for(i in ncol(coldat)){
                    if(all(is.na(coldat[,i])) || 
                       all(coldat[,i]==".") ||
                       all(coldat[,i]=="")){
                        coldat <- coldat[,-i]
                    }
                }
            } 
            #### Check if empty ####
            if(ncol(coldat)==0 | nrow(coldat)==0){
                return(NULL)
            } else {
                ## keeps colnames unchanged
                data.table::as.data.table(coldat)
            } 
        }) ## end lapply
        
        #### Remove NULL entries ####
        nulls <- unname(unlist(lapply(tmp,is.null)))
        tmp <- tmp[!nulls]
        n <- n[!nulls]
        ## Each element can potentially have >1 column 
        ncols <- unlist(lapply(tmp,ncol))
        tmp <- do.call(cbind, tmp)
        #### Add column names ####
        if(!is.null(tmp)){
            if(isTRUE(add_sample_names)){
                colnames(tmp) = paste(rep(n, times = ncols), 
                                      colnames(tmp),sep = "_") 
            } else {
                colnames(tmp) <- rep(n, times = ncols)
            } 
        }
    } 
    df <- cbind(df, tmp) 
    remove(tmp)
    
    ## ----- Post-processing ----- ####
    #### Only keep unique rows ####
    #### Remove duplicated columns ####
    if(isTRUE(unique_cols)){
        drop_duplicate_cols(dt = df)
    } 
    #### Remove any remaining empty columns #####
    if(isTRUE(drop_empty_cols)){
        remove_empty_cols(sumstats_dt = df,
                          sampled_rows = sampled_rows,
                          verbose = verbose)
    }
    #### Unlist columns inplace ####
    if(isTRUE(unique_rows) && isFALSE(unlist_cols)){
        messager("Must set unlist_cols=TRUE to use unique_rows=TRUE.",
                 "Setting unlist_cols=TRUE.",v=verbose)
        unlist_cols <- TRUE
    }
    if(isTRUE(unlist_cols)){ 
        unlist_dt(dt = df,
                  verbose = verbose)
    } 
    #### Remove duplicated rows #### 
    if(isTRUE(unique_rows)){
        df <- drop_duplicate_rows(dt = df, 
                                  verbose = verbose)
    }  
    #### Report ####
    if(verbose) methods::show(round(difftime(Sys.time(),t1),1))
    messager("VCF data.table contains:",
             formatC(nrow(df),big.mark = ","),"rows x",
             formatC(ncol(df),big.mark = ","),"columns.",
             v=verbose)
    return(df)
}
