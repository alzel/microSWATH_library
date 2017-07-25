library(tidyverse)
library(stringr)
library(forcats)

fun_name = "comparing_libraries"
figures_dir = "./figures"

libraries = c("./results/2017-05-05/fractionation_irt_cons_openswath_irt_corrected.csv",
              "./results/2017-05-05/umpire_34_cons_openswath_irt_corrected.csv")
              
proteins <-  unique(library_fractions$ProteinName)

createAbsoluteProteins = function() {
  absolute_data.raw <- readxl::read_xls("./data/2016-09-05/nature02046-s2.xls",sheet = 1) %>% data.frame() 
  absolute_data.raw <- absolute_data.raw[which(grepl(pattern = "[0-9]+\\.?[0-9]+", x = absolute_data.raw$Protein.Molecules.Cell)),] 
  absolute_data.raw$Protein.Molecules.Cell <- as.numeric(as.character(absolute_data.raw$Protein.Molecules.Cell))
  dataset1 <- absolute_data.raw %>% 
    dplyr::select(ORF, Protein.Molecules.Cell) %>% 
    mutate(dataset = "Ghaemmaghami") %>% 
    rename(abundance = Protein.Molecules.Cell)
  
  absolute_data.raw <- readxl::read_xlsx("./data/2016-09-05/nmeth.2834-S2.xlsx", sheet = "S.cerevisiae") %>% data.frame()
  
  absolute_data.raw <- absolute_data.raw[which(grepl(pattern = "[0-9]+\\.?[0-9]+", x = absolute_data.raw$Copy.number)),] 
  absolute_data.raw$Copy.number <- as.numeric(as.character(absolute_data.raw$Copy.number))
  
  dataset2 <- absolute_data.raw %>% 
    dplyr::select(ORF, Copy.number) %>% 
    mutate(dataset = "Kulak") %>% 
    rename(abundance = Copy.number)
  dataset <- bind_rows(dataset1, dataset2)
  absolute_dataset <- droplevels(dataset[dataset$ORF != "",])

  return(absolute_dataset)
}


absolute_dataset <- createAbsoluteProteins()

getProteins <- function(proteins) {

  res <- map(proteins,  .f = function(x) {
    (str_split(string = x, pattern = "/" ) %>% unlist())[-1]
  }) %>% unlist()
  
  return(res )
}
  

library_proteins <- bind_rows(map(.x = libraries, .f = function(x) {
      tmp_data = read_delim(x, delim = "\t")
      tmp_data$file = basename(x)
      
      res <- tibble(file = basename(x),
                    ORF = getProteins(unique(tmp_data$ProteinName)))
      return(unique(res))
    } ))

library_proteins$lib_name = str_replace(string = library_proteins$file, pattern = "(^[a-z]+).*", replacement = "\\1")


library_proteins <- library_proteins %>% filter(!str_detect(ORF, "reverse_")) %>% 
  mutate(present = 1) %>% 
  complete(file, ORF) %>% # makes complete cases to introduce NA that are not present in one of the files 
  mutate(present = ifelse(is.na(present), 0, present)) %>%
  group_by(ORF) %>%
    mutate(inBoth = sum(present))
  

p <- library_proteins %>% left_join(absolute_dataset) %>% filter(!is.na(dataset)) %>% filter(present == 1) %>%
  ggplot(aes(x = log(abundance), fill=lib_name)) +
    geom_density(alpha = 0.25) +
    facet_grid(~dataset, scales = "free") +
    theme(aspect.ratio = 8/5)
    

file_name = paste(fun_name, "png", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p, filename=file_path, width = 11.69, height = 8.27, dpi = 150)

