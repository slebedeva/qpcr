#' An import qpcr function.
#' 
#' This function will import raw exported txt file from BioRad machine (make sure to check all fields when exporting!) and return a table of values including all fields and NA.
#' @param workdir Your working directory where to look for txt files. Defaults to current directory.
#' @export
#' @examples 
#' importqpcr()
 
importqpcr <- function(workdir="."){
  files <- list.files(path = workdir, pattern = ".txt")
  list_of_tabs <- lapply(1:length(files), function(i){
    file <- readLines(files[[i]])
    idx <- which(grepl("Well",file))[2] #Well is mentioned twice, the second one where the table starts
    file <- file[idx:length(file)]
    tc <- textConnection(file)
    processed <- read.table(tc, sep = "\t", dec=".", header=T, blank.lines.skip = TRUE, na.strings = "N/A")
    tab <- processed
    close(tc) 
    tab <- as.data.frame(tab)
    #typeof(tab$Cq) #check data type #double
    if( typeof(tab$Cq)!="double" ) stop('Please check your input!')
    do.call(what = "<-", args = list(paste0("tab", as.character(i)), tab))
  }) #end of listing subfunction
} # end of import function
