
rm(list=ls())
output_dir = "./R/objects"
figures_dir = "./figures"
dir.create(output_dir)
dir.create(figures_dir)

fun_name = "presence_inalldata"

library(plyr);library(dplyr)
library(reshape2)
library(ggplot2)
library(sva)
library(tidyr)

lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

plots.list = list()


load("./R/objects/dataset.YSBN11.raw._load_.RData")

dataset.YSBN11.raw <- tbl_df(dataset.YSBN11.raw)
dataset.YSBN11.raw$E.Name <- as.character(dataset.YSBN11.raw$E.Name)  
dataset.YSBN11.raw$E.Name <- sub(x = dataset.YSBN11.raw$E.Name, pattern = "_(exhaustion|frac)", replacement = "\\1")
# dataset.YSBN11.raw <- dataset.YSBN11.raw %>% 
#   separate(E.Name, into = c("analysis_date", "sample_type", "aquisition_type", "window_size", "library" ), remove = F) %>% 
#   dplyr::select(-analysis_date, -sample_type, -aquisition_type)

dataset.YSBN11 <- dataset.YSBN11.raw %>% 
  dplyr::select(E.Name, R.Label, EG.Label, EG.StrippedSequence, EG.Qvalue) %>% 
  distinct() %>%
  group_by(E.Name, R.Label, EG.Label) %>% 
  filter(EG.Qvalue == min(EG.Qvalue))
       
dataset.YSBN11.f <- dataset.YSBN11 %>% filter(EG.Qvalue < 0.01)


EG.Label.presence <- plyr::ddply(.data = dataset.YSBN11.f, 
                                 .variables = .(E.Name),
                                 .fun = function(x) {
                                      #x = dataset.YSBN11.f %>% filter(E.Name == "151204_YSBN11_mSWATH_34_BiognosysLib")
                                      #x %>% arrange(R.Label, EG.StrippedSequence) %>% View()
                                      
                                      x.tmp <- reshape2::dcast(x, "EG.Label~R.Label", value.var = "EG.Qvalue")
                                      x.tmp.long <- melt(x.tmp, id.vars = "EG.Label")
                                      return(x.tmp.long)
                                    })

EG.Label.presence <- tbl_df(EG.Label.presence %>% mutate(isPresent = ifelse (!is.na(value), T, F)))
View(EG.Label.presence)
EG.Label.presence.stats <- EG.Label.presence %>% group_by(E.Name, EG.Label) %>% 
  summarise(n_samples = sum(isPresent)) %>% 
  group_by(E.Name) %>% View()
  mutate(total = n()) %>% group_by(E.Name, n_samples) %>%
  summarise(sample_count = n(),
            sample_fraction =sample_count/total[1])

EG.Label.presence.stats <- EG.Label.presence.stats %>% 
  separate(E.Name, into = c("analysis_date", "sample_type", "aquisition_type", "window_size", "library" ), remove = F) %>% 
  dplyr::select(-analysis_date, -sample_type, -aquisition_type)
library(ggplot2)
toPlot <- EG.Label.presence.stats

p <- ggplot(toPlot, aes(x = factor(n_samples, levels = rev(unique(n_samples))), y = sample_count, fill = library)) +
       geom_bar(stat = "identity", position = "dodge") +
       facet_wrap(~window_size, scales = "free") +
       ylab("Total peptides quantified, FDR<0.01 ") +
       xlab("Presence in N injections")
plots.list = lappend(plots.list, p)

p <- ggplot(toPlot, aes(x = factor(n_samples, levels = rev(unique(n_samples))), y = sample_fraction, fill = library)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~window_size, scales = "free") +
  ylab("Percentage of total peptides quantified, FDR<0.01 ") +
  xlab("Presence in N injections")
plots.list = lappend(plots.list, p)



dataset.YSBN11.raw %>% filter(E.Name == "151204_YSBN11_mSWATH_29_BiognosysLib", FG.Id == "_LRTDETLR_.2") %>% View()


file_name = paste("supplementery", fun_name, sep = ".")
file_path = paste(figures_dir, file_name, sep="/")

lapply(seq_along(plots.list) , 
       function(x) {
         
         tryCatch({
           p <- plots.list[[x]]
           scale = 1
           if (length(p$toScale) != 0 && p$toScale == T  ){
             scale = 2
           }
           ggplot2::ggsave(filename = paste(file_path, x , "pdf", sep = "."), device = NULL,
                           plot = p, width = 210 , height = 297, units = "mm", scale = scale)
           
           ggplot2::ggsave(filename = paste(file_path, x , "png", sep = "."), device = NULL,
                           plot = p, width = 210 , height = 297, dpi = 150, units = "mm", scale = scale)
           
         }, error = function(e) {
           message(paste("Plot", "x", "sucks!" ))
           return(NULL)
         }, finally = {
           message(paste("processed plot", x))
         })
       })







