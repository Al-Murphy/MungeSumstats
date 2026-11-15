#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if(length(args)!=4) stop('mode panel reference_path output_directory')
mode <- args[[1]];panel_name <- args[[2]];reference_path <- args[[3]];out <- args[[4]]
base <- '/private/tmp/gwas-harvest-stage2'
.libPaths(c(if(mode=='legacy') file.path(base,'Rlib/legacy') else Sys.getenv('MSS_NEW_RLIB',file.path(base,'Rlib')),file.path(base,'Rlib'),'/private/tmp/harvest-Rlib',.libPaths()))
library(MungeSumstats);library(data.table);setDTthreads(1)
assets <- file.path(base,'assets/mss/functionality')
dir.create(out,recursive=TRUE)
input <- file.path(assets,paste0(panel_name,'-input.tsv'))
metadata <- list(mode=mode,panel=panel_name,input=input,source_md5=unname(tools::md5sum(input)))
if(mode %in% c('new','trusted')) {
    columns <- list(chr='CHR',pos='BP',id='SNP',other='A1',effect_allele='A2',effect='BETA',eaf='FRQ')
    if(mode=='trusted') columns <- c(columns,list(source_ref='A1',source_alt='A2'))
    elapsed <- system.time(result <- harmonise_file(input,reference_path,file.path(out,'edition'),frequency_type='EAF',
                         columns=columns,source_strand=if(mode=='trusted') '+' else NULL))
    metadata$source_rows <- result$source_rows
    metadata$retained_rows <- result$normalized_rows
    metadata$status_counts <- result$status_counts
    metadata$implementation_version <- if(is.null(result$implementation_version)) 'mss-preservation-2' else result$implementation_version
    metadata$loaded_dll_md5 <- unname(tools::md5sum(getLoadedDLLs()[['MungeSumstats']][['path']]))
    metadata$source_strand <- result$source_strand
} else {
    # Matched loader adapter, SNPs only, as documented in the scale benchmark.
    # Legacy SNPlocs has no indels; adding indels to it would change legacy QC.
    sites <- fread(file.path(reference_path,'source.tsv'))
    multis <- sites[grepl(',',alt,fixed=TRUE)]
    if(nrow(multis)) {
        expanded <- multis[,list(alt=unlist(strsplit(alt,',',fixed=TRUE))),by=c('chr','pos','ref','id')]
        sites <- rbindlist(list(sites[!grepl(',',alt,fixed=TRUE),list(chr,pos,ref,alt,id)],expanded),use.names=TRUE)
    }
    sites <- sites[nchar(ref)==1 & nchar(alt)==1]
    map <- Biostrings::IUPAC_CODE_MAP
    codes <- setNames(names(map),vapply(map,function(s) paste(sort(strsplit(s,'')[[1]]),collapse=''),character(1)))
    refs <- sites[,list(seqnames=chr[[1]],pos=pos[[1]],ref_allele=ref[[1]],
                       alt_alleles=paste(unique(alt),collapse=','),
                       alleles_as_ambig=unname(codes[paste(sort(unique(c(ref,alt))),collapse='')])),by=id]
    setnames(refs,'id','SNP');setkey(refs,SNP)
    ns <- asNamespace('MungeSumstats')
    loader <- function(snps,...) {x<-refs[unique(snps),nomatch=NULL];setkey(x,SNP);x}
    unlockBinding('load_ref_genome_data',ns);assign('load_ref_genome_data',loader,ns);lockBinding('load_ref_genome_data',ns)
    metadata$source_rows <- nrow(fread(input))
    elapsed <- system.time(result <- format_sumstats(input,ref_genome='GRCh37',dbSNP=144,
                on_ref_genome=FALSE,infer_eff_direction=FALSE,check_dups=FALSE,
                strand_ambig_filter=FALSE,frq_is_maf=FALSE,sort_coordinates=FALSE,
                nThread=1,return_data=TRUE,save_path=file.path(out,'legacy.tsv')))
    metadata$retained_rows <- nrow(result)
    metadata$loader_scope <- 'exact matched SNP references; no indels in legacy loader; reference ID imputation disabled via all valid rsID selection'
}
metadata$timing <- elapsed
saveRDS(metadata,file.path(out,'receipt.rds'))
write.table(data.frame(key=names(metadata),value=vapply(metadata,function(x) paste(x,collapse=';'),character(1))),
            file.path(out,'receipt.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
cat('Functionality complete:',mode,panel_name,metadata$source_rows,metadata$retained_rows,'\n')
