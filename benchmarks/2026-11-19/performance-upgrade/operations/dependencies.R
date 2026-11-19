# Freeze repository sources on a separately admitted image builder, then replay
# the same archives offline. Each attempted install writes its own retained log.
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
mode <- args[[1]]
description <- args[[2]]
root <- args[[3]]
argument <- args[[4]]
sha <- function(path) {
    x <- system2("sha256sum", shQuote(path), stdout=TRUE)
    if (!is.null(attr(x,"status")) || !length(x)) stop("sha256sum failed")
    sub(" .*", "", x[[1]])
}
if (mode == "freeze") {
    if (file.exists(root) || !dir.create(root)) stop("Need a new archive directory")
    repos <- readLines(argument, warn=FALSE)
    if (!length(repos) || any(!grepl("^https://[^/@[:space:]]+/", repos)) ||
        any(grepl("[@?#]", repos))) stop("Supply public repository URLs, one per line")
    db <- utils::available.packages(repos=repos, type="source")
    d <- read.dcf(description)
    parts <- unlist(strsplit(paste(d[,intersect(colnames(d), c("Imports","LinkingTo","Depends","Suggests"))],
                                 collapse=","), ","))
    roots <- setdiff(trimws(sub("\\s*\\(.*$", "", parts)), c("R",""))
    # Priority base packages are furnished by the exact R engine. Recommended
    # packages in the system library are also part of the immutable image.
    system_pkgs <- rownames(utils::installed.packages(lib.loc=.Library))
    roots <- setdiff(roots, system_pkgs)
    missing <- setdiff(roots, rownames(db))
    if (length(missing)) stop("Unavailable package(s): ", paste(missing, collapse=", "))
    edges <- c("Depends","Imports","LinkingTo")
    closure <- unique(c(roots, unlist(tools::package_dependencies(roots, db, which=edges,
                                                                  recursive=TRUE))))
    pending <- sort(setdiff(closure, system_pkgs))
    direct <- tools::package_dependencies(pending, db, which=edges)
    ordered <- character()
    while (length(pending)) {
        ready <- pending[vapply(direct[pending], function(x) !any(x %in% pending), logical(1))]
        if (!length(ready)) stop("Dependency cycle")
        ordered <- c(ordered, ready)
        pending <- setdiff(pending, ready)
    }
    rows <- lapply(ordered, function(pkg) {
        archive <- paste0(pkg,"_",db[pkg,"Version"],".tar.gz")
        url <- paste0(db[pkg,"Repository"],"/",archive)
        path <- file.path(root,archive)
        utils::download.file(url,path,mode="wb",quiet=TRUE)
        data.frame(package=pkg,version=db[pkg,"Version"],archive=archive,source=url,
                   bytes=file.info(path)$size,sha256=sha(path))
    })
    utils::write.table(do.call(rbind,rows),file.path(root,"lock.tsv"),sep="\t",
                       quote=FALSE,row.names=FALSE)
} else if (mode == "install") {
    # Pin lock.tsv externally across transfer before calling this phase.
    lib <- argument
    if (file.exists(lib) || !dir.create(lib)) stop("Need a new shared library")
    logdir <- paste0(lib, "-install-logs")
    if (!dir.create(logdir)) stop("Need a new log directory")
    lock <- read.delim(file.path(root,"lock.tsv"),colClasses="character",check.names=FALSE)
    stopifnot(identical(names(lock),c("package","version","archive","source","bytes","sha256")),
              !anyDuplicated(lock$package), all(grepl("^[A-Za-z][A-Za-z0-9.]+$",lock$package)),
              all(lock$archive == paste0(lock$package,"_",lock$version,".tar.gz")),
              all(basename(lock$archive) == lock$archive), all(grepl("^[a-f0-9]{64}$",lock$sha256)))
    Sys.setenv(R_LIBS_USER=lib,R_LIBS_SITE=lib,MAKEFLAGS="-j1")
    for (i in seq_len(nrow(lock))) {
        p <- file.path(root,lock$archive[i])
        if (!identical(sha(p),lock$sha256[i])) stop("Dependency archive differs: ", p)
        log <- file.path(logdir,paste0(lock$package[i],".log"))
        rc <- system2(file.path(R.home("bin"),"R"),
                      c("CMD","INSTALL","--preclean",shQuote(paste0("--library=",lib)),shQuote(p)),
                      stdout=log,stderr=log)
        if (rc != 0L) stop("Dependency install failed: ",lock$package[i])
        got <- utils::packageDescription(lock$package[i],lib.loc=lib,fields="Version")
        stopifnot(identical(as.character(got),lock$version[i]))
    }
} else stop("Expected freeze or install")
