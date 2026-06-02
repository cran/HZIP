## ## for creating a package
.package.Name <- "HZIP"

##.First.lib <- function(lib,pkg){
##  library.dynam(.package.Name,
##                pkg,
##                lib)
##}

##.Last.lib <- function(libpath){
##  library.dynam.unload(chname="HZIP",libpath=libpath)
##}

.onAttach <- function(...) {
  # Texto base
  msg_lines <- c(
    "Classes and Methods for R originally developed in the",
    "Complex Statistical Modeling Laboratory (CoSMo)",
    "Department of Statistics",
    "Federal University of Bahia, Brazil (2025),",
    "by and under the direction of",
    "Jalmar M. F. Carrasco and Lizandra C. Fabio,",
    "with contributions from collaborators and students.",
    "Main functions: hzip, residuals, envelope."
  )

  # Calcula largura máxima e define margem
  width <- max(nchar(msg_lines))
  border <- strrep("*", width + 6)  # 3 espaços de cada lado

  # Monta texto com moldura
  boxed <- c(
    border,
    paste0("*  ", format(msg_lines, width = width), "  *"),
    border
  )

  # Exibe com packageStartupMessage
  packageStartupMessage(paste(boxed, collapse = "\n"))
}

.onUnload <- function(libpath){
  library.dynam.unload("HZIP",libpath=libpath)
}
