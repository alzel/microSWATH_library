#!/usr/bin/env Rscript
# analysing DDA data

library(tidyverse)
library(stringr)
library(forcats)


fun_name = "20170215_dda_analysis"
figures_dir = "./figures"
dir.create(figures_dir, recursive = T)

data.raw <- read_csv("./results/2017-02-15/iprophet_psm_protFDR0.01_t_1.07.csv")
data.raw <- data.raw %>% 
  separate(col = scan, into = c("file", "scan"), sep = "\\.", extra = "merge") %>% 
  mutate(file = str_replace(file, "lowmicro_", "lowmicro"),
         file = str_replace(file, "-.*?$", "")) 
write_tsv(x = data.raw, path = "./results/2017-02-15/iprophet_psm_protFDR0.01_t_1.07_tidy.csv")

p.counts <- data.raw %>% group_by(file) %>% 
  filter(mFDR < 0.01) %>%
  summarize(n_peptides = n_distinct(pep),
            n_proteins = n_distinct(prot)) %>%
  gather(n_peptides, n_proteins, key = "variable", value = "count") %>%
  ggplot(aes(x = fct_reorder(file, count), y = count, fill = variable ) )+
        geom_bar(stat = "identity", position = "dodge") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

file_name = paste(fun_name, "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p.counts, filename=file_path, width = 11.69, height = 8.27)


