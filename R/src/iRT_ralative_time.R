#!/usr/bin/env Rscript
# getting relative iRTs for microSWATH library


rm(list=ls())
plots.list = list()

fun_name = "iRT_relative_time"
figures_dir = "./figures"

dir.create(figures_dir, recursive = T)

library(dplyr)
library(ggplot2)

v4lib_QC_apex <- read.delim("./data/2015-12-09/v4lib_QC_apex.tsv")

iRT_fasta = paste0(unique(v4lib_QC_apex$Peptide.Sequence),collapse = "")
rev_iRT_fasta = paste(rev(unlist(strsplit(iRT_fasta, split = ""))), collapse = "")

thr=2

normalized = v4lib_QC_apex %>% group_by(File.Name) %>% mutate(n = length(apex)) %>% filter(n == 11) %>% mutate(iRT = (apex - sort(apex)[2])/(max(apex) - sort(apex)[2]))
normalized = normalized %>% group_by(Peptide.Sequence) %>% mutate( Z_apex = (apex - mean(apex))/sd(apex))
normalized$toRemove = ifelse(abs(normalized$Z_apex) > thr, 1, 0)


p = ggplot(normalized, aes(x = File.Name, y=iRT, colour=Peptide.Sequence)) + 
      geom_point() +
      geom_point(data = normalized%>%filter(toRemove == 1), shape=4, colour="red", size=5) +
      theme(axis.text.x = element_text(angle=90, vjust=1),
            aspect.ratio = 5/8)

file_name = paste(fun_name, "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p, filename=file_path, width = 11.69, height = 8.27)


iRT.txt = normalized %>% filter(toRemove != 1 ) %>% group_by(Peptide.Sequence) %>% summarise(mean = round(mean(iRT, na.rm=T)*100,4), sd = sd(iRT)) %>% arrange(mean) %>% dplyr::select(Peptide.Sequence, mean)

file_path = "./results/2015-12-05/iRT_microSWATH.txt"
write.table(x = iRT.txt, file = file_path, row.names = F, col.names = F, quote = F, sep = "\t")



#fixing irt in the libraries
irt.dataset  = normalized %>% filter(toRemove != 1 ) %>% group_by(Peptide.Sequence) %>% 
               summarise(meanIRT = round(mean(iRT, na.rm=T)*100,4), sd = sd(iRT), meanRT = mean(apex, na.rm=T)*60) %>% 
               arrange(meanIRT) %>% dplyr::select(Peptide.Sequence, meanIRT, meanRT)



input_path = "./results/2016-01-04/spectrast_results"
output_dir = "./results/2016-01-04/final_libraries"
dir.create(output_dir)


filesToProcess = dir(path=input_path, pattern = "*openswath.csv", recursive=F, full.names = T)

for (i in 1:length(filesToProcess)) {
  input_file = filesToProcess[i]
  lib.raw <- read.delim(input_file)
  tmp.dataset = data.frame(rt = lib.raw$Tr_recalibrated)
  tmp = cbind(irt.dataset, data.frame (rt = lib.raw[match(irt.dataset$Peptide.Sequence, lib.raw$PeptideSequence), "Tr_recalibrated"]))
  lm.fit = lm(meanIRT~rt, data = na.omit(tmp))
  
  lib.raw$Tr_recalibrated = predict(lm.fit, tmp.dataset)
  output_file = sub(pattern = ".csv", x = basename(input_file), replacement = "_iRT.csv")
  tmp.file_path = paste(output_dir, output_file, sep = "/")
  write.table(x = lib.raw, file = tmp.file_path, row.names = F, col.names = T, quote = F, sep = "\t")
}




