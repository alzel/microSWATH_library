source("./R/boot.R")
rm(list=ls())
output_dir = "./R/objects"
dir.create(output_dir)
suffix = "_load_"

##################################################
###loading proteome data  
##################################################


## loads consequtive injections, processed with different spectranaut libraries
load_precision_data = function() {
  
  input_path = "./data/2016-09-22/microSWATH/benchmark/analysis"
  message("Loading proteome")
  
  filesToProcess = dir(path=input_path, pattern = "Report_Standard_SWATH_Fragments_v2\\ \\(Normal\\).xls", recursive=T)
  
  #filesToProcess = filesToProcess[grep(pattern="*.xls", x=filesToProcess)]
  #pattern.p = "batch([0-9]+)*.*"
  #pattern.p = "KL_Try_Spectronaut_Batch02_([0-9]+)_([A-Za-z0-9 ]+)_([0-9]+)_([0-9]+).xls"
  
  matches = stringr::str_match_all(pattern=".*YSBN11.*", filesToProcess)
  
  read_peptides = function(x) {
    file_name = paste(input_path,x[[1]], sep="/") 
    table = read.table(file_name, sep="\t", header=T)
    #table$acquisition_batch = factor(rep(as.integer(x[[2]]), nrow(table))) #Spectronaut batch
    #table$type  = factor(rep(x[[3]], nrow(table)))
    #table$Sp_date =  factor(rep(x[[4]], nrow(table))) #Spectranaut date
    #table$Sp_time =  factor(rep(x[[5]], nrow(table))) 
    table$file =  as.character(x[[1]])
    return(table)
  }
  
  
  file.list = lapply(matches, FUN=read_peptides)
  dataset.YSBN11.raw = tbl_df(do.call(rbind.data.frame, file.list))
  
  file_name = paste("dataset.YSBN11.raw", suffix, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(dataset.YSBN11.raw,file=file_path)

}






