rm(list = ls())
library(tidyverse)
library(forcats)
library(stringr)

#analyzing proteomes obtaned from Y3K paper
fun_name = "coons_stats"
figures_dir = "./figures"
dir.create(figures_dir, recursive = T)

input_path = "./data/2017-05-19/"
filesToProcess = dir(path=input_path, pattern = "^Protein_.*?txt", recursive=T)
filesToProcess = grep(pattern="", filesToProcess, value=T)

pattern.p = "Protein_([A-Za-z]+)_[0-9]+.txt$"

matches = stringr::str_match_all(pattern=pattern.p, filesToProcess)

read_data = function(x) {
  
  file_name = paste(input_path,"/Protein_", x[[2]],"/", x[[1]], sep="") 
  #file_name = './data/2017-05-19//Protein_Respiration/Protein_Respiration_9.txt'
  tmp_data = read_delim(file_name, delim = "\t", trim_ws = F, progress = T)
  tmp_data <- tmp_data %>% dplyr::select(-matches("X\\d+")) 
  
  tmp_data <- tmp_data %>% dplyr::select(`Protein IDs`, matches("LFQ intensity")) %>% gather(mutant, "LFQ", -`Protein IDs`)
  tmp_data$condition <- x[[2]]
  tmp_data$file <- x[[1]]
  
  
  return(tmp_data)
}

data.list = lapply(matches, FUN=read_data)
dataset <- bind_rows(data.list)
dataset <- dataset %>% mutate(mutant = str_replace(mutant, pattern = "LFQ intensity ", replacement = "" ),
                              mutant = paste(mutant, file, sep = "|"))
dataset <- dataset %>% rename(protein_id = `Protein IDs`)


toPlot <- dataset %>% 
  group_by(mutant, condition, file) %>% 
  summarise(total_signal = sum(LFQ, na.rm = T)) %>%
  ungroup() %>%
  arrange(file) %>%
  mutate(mutant_id  = seq_along(mutant))

p = toPlot %>%
  ggplot(aes(y = total_signal, x = mutant_id)) +
    geom_point(aes(col = file, shape = condition)) +
    geom_point(data = toPlot %>% filter(grepl("WT", mutant)), aes(col = "red", size = 4 , shape=condition))

file_name = paste("batches", fun_name, "png", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p, filename=file_path, width = 11.69, height = 8.27, dpi = 150)


toPlot <- dataset %>% filter(grepl("WT", mutant), !is.na(LFQ)) %>% 
  group_by(protein_id, file) %>% 
  mutate(total_rep = n()) %>% 
  filter(total_rep >=3 ) %>% 
  summarise(CV = sd(LFQ)/mean(LFQ))

p <- ggplot(toPlot, aes(x = CV)) +
    geom_density() +
    geom_vline(xintercept = median(toPlot$CV)) +
    geom_text(x = 0.2, y = 2, label= median(toPlot$CV)) +
    xlim(0,1)

file_name = paste("CVs", fun_name, "png", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p, filename=file_path, width = 11.69, height = 8.27, dpi = 150)








