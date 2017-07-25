rm(list = ls())
library(tidyverse)
library(forcats)
library(stringr)

fun_name = "peak_shape"
figures_dir = "./figures"


lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

plots.list = list()

read_delim("./results/2017-03-18/20170410_104253_nselevse_Report.xls", delim = "\t", col_types = cols(E.Name = col_character()))


data_raw <- bind_rows(read_delim("./results/2017-03-18/20170410_104253_nselevse_Report.xls", delim = "\t", col_types = cols(E.Name = col_character())),
                      read_delim("./results/2017-03-18/20170407_140301_151204_Report.xls", delim = "\t", col_types = cols(E.Name = col_character())),
                      read_delim("./results/2017-03-18/20170318_150344_151204_Report.xls", delim = "\t", col_types = cols(E.Name = col_character())))
                      


x.df = data_raw %>% filter(R.FileName == "nselevse_L120412_001_SW")
y.df = data_raw %>% filter(R.FileName == "151204_YSBN11_mSWATH_29_01")


common_precursors <- (semi_join(x.df, y.df, by = c("EG.PrecursorId")) %>% dplyr::select(EG.PrecursorId))$EG.PrecursorId


toPlot <- data_raw %>% filter(EG.PrecursorId %in% common_precursors) %>% 
  mutate(ratio.a_h = F.PeakArea/F.PeakHeight) 
  
p <- ggplot(toPlot, aes(y = ratio.a_h, x = R.FileName)) +
  geom_boxplot() +
  scale_y_continuous(limits = quantile(toPlot$ratio.a_h, c(0.1, 0.95))) +
  scale_x_discrete( labels = c("micro/29x16", "micro/34x25", "nano/34x25")) +
  xlab("Chromatography/SWATH regime") +
  ylab("Peak Area to Peak height ratio")



file_name = paste(fun_name, "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p, filename=file_path, width = 11.69, height = 8.27)

file_name = paste(fun_name, "png", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(plot = p, filename=file_path, width = 11.69, height = 8.27, dpi = 150)
