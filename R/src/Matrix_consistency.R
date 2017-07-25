rm(list = ls())
library(tidyverse)
library(forcats)
library(stringr)

fun_name = "Matrix_consistenscy"
figures_dir = "./figures"

lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

plots.list = list()


READDATA <- function(sheet) {
  data.raw <- readxl::read_excel("./data/2017-03-20/nbt.3683-S4.xlsx", sheet = sheet)
  names(data.raw) <- str_replace(string = names(data.raw), pattern = "∆", replacement = "")
  names(data.raw) <- str_replace(string = names(data.raw), pattern = " ", replacement = "_")
  
  data.f <- data.raw %>% filter(Molecule_Type == "Protein") %>% 
    gather(knockout, LFQ,  -Molecule_Name, -Molecule_Type) %>% 
    mutate(knockout = paste(knockout, sheet, sep = "|"))
}



makeNA <- function(data.f, fraction) {
  #setts fraction of good LFQ to NA data.f
  data.f.present <- data.f %>% filter(!is.na(LFQ))
  data.f.absent <- data.f %>% filter(is.na(LFQ))
  data.f.unimputed <- data.f.present %>% sample_frac(size = 1-fraction)
  data.f.unimputed.NA <- anti_join(data.f.present, data.f.unimputed, by =c("Molecule_Name", "Molecule_Type", "knockout" )) %>% mutate(LFQ = NA)
  data.NA <- bind_rows(list(data.f.unimputed, data.f.unimputed.NA, data.f.absent))
  return(data.NA)
}


getSliced <- function(data.f) {
  ko_list <- sample(unique(data.f$knockout))
  get_data <- function(i, ko_list) {
    ko_list[i]
  }
  
  my_filter <- function(x, data.f) {
    data.f %>% filter(knockout %in% x)
  }
  
  indices <-  map(1:length(ko_list), seq)
  
  names_idx <- map(indices, get_data, ko_list)
  
  data.sliced <- bind_rows(map(names_idx, my_filter, data.f), .id = "id")
  data.sliced.s <- data.sliced %>% group_by(id, Molecule_Name) %>% summarise(present = !any(is.na(LFQ)))
  return(data.sliced.s)  
}


data1 <- READDATA(1)
data4 <- READDATA(4)
data <- bind_rows(list(data1, data4))

#removing fraction (unimputing imputed) as reported in http://www.nature.com/nbt/journal/v34/n11/full/nbt.3683.html 
set.seed(123)
data1.sim <- makeNA(READDATA(1), fraction =  4.05/100)
data4.sim <- makeNA(READDATA(4), fraction = 4.53/100)
data.sim <- bind_rows(list(data1.sim, data4.sim))

set.seed(123)
data.sliced <- getSliced(data)
data.sliced$type = "imputed"

set.seed(123)
data.sim.sliced <- getSliced(data.sim)
data.sim.sliced$type = "unimputed"

toPlot <- bind_rows(data.sliced, data.sim.sliced) %>% group_by(id, type) %>% summarise(n_present = sum(present))

p <- ggplot(toPlot, aes(x = as.numeric(id), y = n_present)) + 
  geom_line(aes(colour = type)) + 
  xlab("sample") + 
  ylab("Nr of proteins with at least 1 peptide with < FDR 0.01 ") +
  theme(aspect.ratio = 5/8)

file_name = paste(fun_name, "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p, filename=file_path, width = 11.69, height = 8.27)

file_name = paste(fun_name, "png", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p, filename=file_path, width = 11.69, height = 8.27, dpi = 150)


p <- ggplot(toPlot %>% filter(type != "imputed") , aes(x = as.numeric(id), y = n_present)) + 
  geom_line() + 
  xlab("sample") + 
  ylab("Nr of proteins with at least 1 peptide with < FDR 0.01 ") +
  theme_bw(base_family = "Helvetica") +
  theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size=16)) +
  theme(aspect.ratio = 5/8)

file_name = paste(fun_name, 2,"png", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p, filename=file_path, width = 11.69, height = 8.27, dpi = 150)

