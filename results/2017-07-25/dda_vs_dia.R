library(tidyverse)
library(stringr)
rm(list = ls())

output_dir = "./results/2017-07-25/"
figures_dir = "./figures"
fun_name = "dda_vs_dia"
dir.create(output_dir)
dir.create(figures_dir)


loadData = function() {
  
  
  data_dir = "./data/2016-01-08/"
  paths <- dir(data_dir, pattern = ".*peptides\\.tsv$", full.names = TRUE)
  names(paths) <- basename(paths)
  dataset.raw = plyr::ldply(paths, read.delim, stringsAsFactors = T)
 
  dataset.raw$file = dataset.raw$.id
  dataset.raw$batch = as.factor(sub(pattern = ".*batch_(\\d+).*", replacement = "\\1" , x = dataset.raw$.id))
  
  #data.raw.f = data.raw %>% dplyr::select(ID, uniProt_ID, StrippedSequence, Sample, NormArea, batch)
  data.raw.f = tbl_df(dataset.raw) %>% dplyr::select(ID, uniProt_ID, StrippedSequence, Sample, NormArea,Qvalue, batch)
 
 
  
  #suffix = paste("_", match.call()[[1]], "_", sep = "")
  
  file_name = paste("data.raw.f", fun_name, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(data.raw.f,file=file_path)
  
  file_name = paste("dataset.raw", fun_name, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(dataset.raw,file=file_path)  
}
#laodData()

## -- preparing DIA data ---- 
load(paste(output_dir, "dataset.raw.dda_vs_dia.RData", sep = "/"))

dataset_dia <- dataset.raw %>% as.tibble() %>% ungroup() %>%
  dplyr::select(uniProt_ID, StrippedSequence, Sample, NormArea, Qvalue, batch) %>% 
  mutate(Sample = paste(Sample, batch, sep = "_")) %>%
  mutate(Qvalue = ifelse(Qvalue > 0.01, NA, Qvalue)) %>%
  filter(!is.na(Qvalue)) %>% arrange(Sample, uniProt_ID) %>% 
  group_by(uniProt_ID, Sample) %>% 
    top_n(1, wt = NormArea) %>% 
  ungroup() %>%
  dplyr::select(uniProt_ID, Sample, NormArea) %>%
  complete(uniProt_ID, Sample)
  

dataset_dia <- dataset_dia %>% 
  rename(protein_id = uniProt_ID,
         mutant = Sample,
         LFQ = NormArea) %>%
  mutate(protein_id = as.character(protein_id),
         mutant = as.character(mutant))
         


length(unique(dataset_dia$mutant))


getSliced <- function(data.f, n) {
  # data.f <- dataset_dia
  # n = 10
  ko_list <- sample(unique(data.f$mutant), n)
  
  # ko_list <- (data.f %>% filter(mutant %in% ko_list) %>% 
  #               group_by(mutant) %>% 
  #               summarise(total = sum(!is.na(LFQ))) %>% 
  #               ungroup() %>% 
  #               arrange(desc(total)) %>% 
  #               dplyr::select(mutant))$mutant
  
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


dataset_dia.sliced <- getSliced(data.f = dataset_dia, n = 331)



### --- clearning DDA data --- 
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

dataset_dda <- dataset %>% ungroup %>% complete(protein_id, mutant)

dataset_dda.sliced <- getSliced(data.f = dataset_dda, n = 1000)

dataset_dda.sliced$type = "DDA"
dataset_dia.sliced$type = "DIA"


dataset_combined <- bind_rows(dataset_dda.sliced, dataset_dia.sliced)


iround <- function(x, interval){
  ## Round numbers to desired interval
  ##
  ## Args:
  ##   x:        Numeric vector to be rounded
  ##   interval: The interval the values should be rounded towards.
  ## Retunrs:
  ##   a numeric vector with x rounded to the desired interval.
  ##
  interval[ifelse(x < min(interval), 1, findInterval(x, interval))]
}



dataset_combined <- dataset_combined %>% 
  mutate(id = as.numeric(id)) %>% 
  group_by(type) %>% 
  mutate(fraction = n_present/max(n_present),
         n_norm = id/max(id))


dataset_combined = dataset_combined %>% mutate(n_norm_rounded = iround(n_norm, c(0.05, 0.1, 0.2, 0.5, 0.9)))

p_fraction <- ggplot(data=dataset_combined, aes(x = n_norm, y = fraction)) + 
  geom_line(aes(colour = type)) +
  geom_point(data = dataset_combined %>% group_by(type, n_norm_rounded) %>% filter(n_norm > 0.05) %>% top_n(1, wt = -n_norm), aes(x = n_norm, colour = type, y=fraction)) +
  ggrepel::geom_text_repel(data = dataset_combined %>% group_by(type, n_norm_rounded) %>% filter(n_norm > 0.05) %>% top_n(1, wt = -n_norm), aes(x = n_norm, y=fraction, colour = type, label = n_present)) +
  xlab("Present in X% samples") + 
  ylab("Fraction of proteins with at least 1 peptide with < FDR 0.01 ") +
  theme(aspect.ratio = 5/8)



p_300 <- ggplot(data=dataset_combined %>% filter(id < 330), aes(x = id, y = fraction)) + 
  geom_line(aes(colour = type)) +
  geom_point(data = dataset_combined %>% filter(id %in% c(1, 10, 20, 50, 100, 200, 300)), aes(x = id, colour = type, y=fraction)) +
  ggrepel::geom_text_repel(data = dataset_combined %>% filter(id %in% c(1, 10, 20, 50, 100, 200, 300)), aes(x = id, y=fraction, colour = type, label = n_present)) +
  xlab("Present in n samples") + 
  ylab("Fraction of proteins with at least 1 peptide with < FDR 0.01 ") +
  theme(aspect.ratio = 5/8)

length(intersect(dataset_dda$protein_id, dataset_dia$protein_id))/length(unique(dataset_dia$protein_id))
