rm(list = ls())
library(tidyverse)
library(forcats)
library(stringr)

#analyzing proteomes obtaned from Y3K paper
fun_name = "Matrix_consistenscy2"
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
  tmp_data = read_delim(file_name, delim = "\t", trim_ws = FALSE, progress = T)
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

dataset<-dataset %>% ungroup %>% complete(protein_id, mutant)

getSliced <- function(data.f, n) {
  #data.f <- dataset
  #n = 100
  ko_list <- sample(unique(data.f$mutant), n)
  
  ko_list <- (data.f %>% filter(mutant %in% ko_list) %>% group_by(mutant) %>% 
    summarise(total = sum(!is.na(LFQ))) %>% 
    ungroup() %>% arrange(desc(total)) %>% dplyr::select(mutant))$mutant
  
  get_data <- function(i, ko_list) {
    ko_list[i]
  }
  
  my_filter <- function(x, data.f) {
    data.f %>% filter(mutant %in% x) %>% group_by(protein_id) %>% 
      summarise(present = !any(is.na(LFQ))) %>% 
      ungroup() %>%
      summarise(n_present = sum(present))
  }
  
  
  indices <-  map(1:length(ko_list), seq)

  names_idx <- map(indices, get_data, ko_list)
  data.sliced <- bind_rows(map(names_idx, my_filter, data.f), .id = "id")
  
  
  return(data.sliced)  
}


dataset.sliced <- getSliced(data.f = dataset, n = 1000)


p <- ggplot(data=dataset.sliced, aes(x = as.numeric(id), y = n_present)) + 
  geom_line() + 
  xlab("sample") + 
  ylab("Nr of proteins with at least 1 peptide with < FDR 0.01 ") +
  theme(aspect.ratio = 5/8)

file_name = paste(fun_name, "png", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p, filename=file_path, width = 11.69, height = 8.27, dpi = 150)

