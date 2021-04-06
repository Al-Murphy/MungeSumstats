#' Downloads and loads the SNP_LOC_DATA file that I stashed at figshare
#'
#' @return SNP_LOC_DATA Has four columns with SNP CHR BP and Build
#'
#' @examples
#' SNP_LOC_DATA <- load_snp_loc_data()
#'
#' @export
#' @importFrom utils download.file
load_snp_loc_data <- function(){
  #TODO IF SUBMITTING TO BIOCONDUCTOR, NEED TO ADD THIS TO EXPERIMENT HUB
  #print("There is no SNP column found within the data. It must be inferred from CHR and BP information.")
  #print("Note: this requires downloading a 300mb file from figshare into a temporary directory")
  #print("the file which is downloaded is created by the build_snp_location_tables function included with this package")
  SNP_LOC_DATA=NA # Because SNP_LOC_DATA is loaded from a file, we need to trick devtools::check() into passing with this
  filePath = sprintf("%s/hssGWAS/data/SNP_LOC_DATA.rda",.libPaths()[1])
  if(!file.exists(filePath)){
    default <- getOption('timeout')
    options(timeout=200)
    start.time <- Sys.time()
    utils::download.file("https://ndownloader.figshare.com/files/21768105",destfile=filePath)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(sprintf("File downloaded in %.0f seconds",time.taken))
    options(timeout=default)
  }
  load(filePath)
  return(SNP_LOC_DATA)
}