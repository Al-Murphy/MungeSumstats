library(Rcpp)
base <- "/private/tmp/gwas-harvest-stage2/results/reference-index-capacity-v1"
out <- commandArgs(trailingOnly = TRUE)[[1L]]
started <- proc.time()
ref <- MungeSumstats::prepare_reference(
    fasta = file.path(base, "reference.fa"),
    variants = file.path(base, "prefix-sites.tsv"),
    path = file.path(base, out), build = "GRCh38",
    version = "capacity-only-incomplete-prefix-no-history-1000000-v1", batch_size = 100000L
)
dll <- getLoadedDLLs()[["MungeSumstats"]][["path"]]
saveRDS(list(elapsed = proc.time() - started, reference = ref,
             dll = dll, dll_md5 = tools::md5sum(dll)),
        file.path(base, paste0(out, "-result.rds")))
print(proc.time() - started)
