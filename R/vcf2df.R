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
#' @param verbose Print messages.
#' 
#' @source https://gist.github.com/zhujack/849b75f5a8305edaeca1001dfb9c3fe9
#' @source 
#' \code{
#' #### VariantAnnotation ####
#' vcf_file <- system.file("extdata", "ALSvcf.vcf",
#'                         package = "MungeSumstats")
#' vcf <- VariantAnnotation::readVcf(file = vcf_file)
#' vcf_df <- MungeSumstats::vcf2df(vcf = vcf)
#' } 
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
#' @keywords internal
#' @importFrom utils type.convert 
#' @importFrom data.table as.data.table
#' @importFrom methods slotNames is
vcf2df <- function(vcf, 
                   add_sample_names=TRUE,
                   add_rowranges=TRUE,
                   verbose=TRUE) {
    requireNamespace("VariantAnnotation")
    requireNamespace("MatrixGenerics")
  
    messager("Converting VCF to data.table.",v=verbose) 
    #### Expand VCF ####
    if (methods::is(vcf,"CollapsedVCF")) {
        message('Expanding VCF first, so number of rows may increase',
                v=verbose)
        vcf <- VariantAnnotation::expand(x = vcf)
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
        yy = data.frame(lapply(yy,type.convert))
        colnames(yy) = paste("ANN",anncols,sep="_")
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
            for(i in ncol(coldat)){
                if(all(is.na(coldat[,i])) || 
                   all(coldat[,i]==".") ||
                   all(coldat[,i]=="")){
                    coldat <- coldat[,-i]
                }
            }
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
        if(isTRUE(add_sample_names)){
            colnames(tmp) = paste(rep(n, times = ncols), 
                                  colnames(tmp),sep = "_") 
        } else {
            colnames(tmp) <- rep(n, times = ncols)
        } 
    } 
    df <- cbind(df, tmp) 
    if(verbose) methods::show(round(difftime(Sys.time(),t1),1))
    return(df)
}
