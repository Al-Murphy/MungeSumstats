# ----- R Package Dockerfile -----
#
# This Dockerfile is designed for developers of any R package stored on GitHub.
#
# It runs several steps:
#   1. Pulls the official bioconductor Docker container (which includes Rstudio).
#   2. Runs CRAN checks on the R package.
#   3. Installs the R package and all of its dependencies (including Depends, Imports, and Suggests).
#
# This Dockerfile should be used with the [dockerhub.yml](https://github.com/neurogenomics/orthogene/blob/main/.github/workflows/dockerhub.yml) workflow file,
# as you must first checkout the R package from GitHub,
# along with several other GitHub Actions.
#
# If the R package passes all checks, the dockerhub.yml workflow will subsequently
# push the Docker container to DockerHub (using the username and token credentials
# stored as GitHub Secrets).
#
# You can then create an image of the Docker container in any command line:
#   docker pull <DockerHub_repo_name>/<package_name>
# Once the image has been created, you can launch it with:
#   docker run -d -e ROOT=true -e PASSWORD=bioc -v ~/Desktop:/Desktop -v /Volumes:/Volumes --rm -p 8788:8787 <DockerHub_repo_name>/<package_name>
# Finally, launch the containerised Rstudio by entering the following URL in any web browser:
#   http://localhost:8788/
#
# The username will be "rstudio" by default,
# and you can set the password to whatever you like,
#
# This DockerFile was partly adapted from the [scFlow Dockerfile](https://github.com/combiz/scFlow/blob/master/Dockerfile).
FROM bioconductor/bioconductor_docker:devel
RUN apt-get update && \
    apt-get install -y \
    git-core \
    libcurl4-openssl-dev \
    libgit2-dev \
    libicu-dev \
    libssl-dev \
    make pandoc \
    pandoc-citeproc \
    zlib1g-dev \
	xfonts-100dpi \
	xfonts-75dpi \
	biber \
	libsbml5-dev \
	qpdf \
	cmake \
	&& apt-get clean \
    && rm -rf /var/lib/apt/lists/*
# Create a buildzone folder named after the R package
# BiocCheck requires the buildzone to have the same name as the R package
ARG PKG
RUN echo $PKG
RUN mkdir -p /$PKG
ADD . /$PKG
WORKDIR /$PKG
# Install dependencies with AnVil (faster)
RUN Rscript -e 'options(download.file.method= "libcurl"); \
                if(!require("BiocManager")) install.packages("BiocManager"); \
                if(!require("AnVIL"))  {BiocManager::install("AnVIL", ask = FALSE)}; \
                AnVIL::install(c("remotes","devtools")); \
                try({remotes::install_github("bergant/rapiclient")}); \
                bioc_ver <- BiocManager::version(); \
                options(repos = c(AnVIL::repositories(),\
                                  AnVIL = file.path("https://bioconductordocker.blob.core.windows.net/packages",bioc_ver,"bioc"),\
                                  CRAN = "https://cran.rstudio.com/"),\
                                  download.file.method = "libcurl", Ncpus = 2); \
                deps <- remotes::dev_package_deps(dependencies = TRUE)$package; \
                AnVIL::install(pkgs = deps,  ask = FALSE); \
                deps_left <- deps[!deps %in% rownames(installed.packages())]; \
                if(length(deps_left)>0) devtools::install_dev_deps(dependencies = TRUE, upgrade = "never");'
# Run R CMD check - will fail with any errors or warnings
# Run Rscript -e 'devtools::check()'
# Run Bioconductor's BiocCheck (optional)
#ARG BIOC
#RUN if [ "$BIOC" = "true" ]; then \
#        Rscript -e 'if(!require("BiocCheck")) AnVIL::install("BiocCheck");\
#                    BiocCheck::BiocCheck(`quit-with-status` = TRUE,\
#                    `no-check-R-ver` = TRUE,\
#                    `no-check-bioc-help` = TRUE);'\
#    fi
# Install R package from source
RUN R -e 'remotes::install_local(upgrade="never")'
RUN rm -rf /$PKG
