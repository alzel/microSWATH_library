library(plyr);library(dplyr)
library(reshape2)
library(ggplot2)
library(sva)
library(tidyr)
rm(list = ls())
output_dir = "./R/objects"
figures_dir = "./figures"
fun_name = "batch_effects"
dir.create(output_dir)
dir.create(figures_dir)



clean = function() {
  
  
  data_dir = "./data/2016-01-08/"
  paths <- dir(data_dir, pattern = ".*peptides\\.tsv$", full.names = TRUE)
  names(paths) <- basename(paths)
  dataset.raw = ldply(paths, read.delim, stringsAsFactors = T)
  
  dataset.raw$.id
  dataset.raw$file = dataset.raw$.id
  dataset.raw$file = dataset.raw$.id
  dataset.raw$batch = as.factor(sub(pattern = ".*batch_(\\d+).*", replacement = "\\1" , x = dataset.raw$.id))
  
  #data.raw.f = data.raw %>% dplyr::select(ID, uniProt_ID, StrippedSequence, Sample, NormArea, batch)
  data.raw.f = tbl_df(dataset.raw) %>% dplyr::select(ID, uniProt_ID, StrippedSequence, Sample, NormArea, batch)
  
  
  
  #suffix = paste("_", match.call()[[1]], "_", sep = "")
  
  file_name = paste("data.raw.f", fun_name, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(data.raw.f,file=file_path)
  
  file_name = paste("dataset.raw", fun_name, "RData", sep=".")
  file_path = paste(output_dir, file_name, sep="/")
  save(dataset.raw,file=file_path)  
}

clean()

MM_dates = read.delim("./data/2016-01-08/MM_dates.tsv")
MM_sample_code = read.delim("./data/2016-01-08/MM_sample_code.tsv")
# MM_sample_code %>% 
#   filter(Batch == 2 | Batch ==3) %>% 
#   group_by(Sample, Batch) %>% 
#   mutate(n = n()) %>%
#   ungroup() %>%
#   arrange(Sample)
# 
# 
# MM_sample_code[grep("mix", ignore.case = T, MM_sample_code$File.Name),]
      

load("./R/objects/data.raw.f.batch_effects.RData")

# all_ids = data.raw.f %>% select(ID) %>% unique()
# 
# data.matrix.long = data.matrix %>% 
#   gather(Sample, NormArea, -uniProt_ID, -ID) 
# 
# selected = sample(all_ids$ID, 50)
# 
# data.matrix.long.selected = data.matrix.long[data.matrix.long$ID %in% selected,]
# 
# data.matrix.long.selected2 = data.matrix.long.selected %>%
#   separate(Sample, c("Sample", "batch"), -3, remove = F) %>%
#   mutate(batch = gsub("_", "", batch)) %>%
#   arrange(uniProt_ID, ID)
# data.matrix.long.selected2$is.Present = ifelse(is.na(data.matrix.long.selected2$NormArea), 0, 1)
# 
# data.matrix.long.selected2 = data.matrix.long.selected2 %>% ungroup() %>% arrange(batch)
# data.matrix.long.selected2$Sample = factor(data.matrix.long.selected2$Sample, levels = unique(data.matrix.long.selected2$Sample))
# p = ggplot(data.matrix.long.selected2, aes(x=Sample, y=ID, fill=factor(is.Present))) +
#     geom_raster()

#MM_dates$Sample = sub(pattern = "\\d+_AK[_]?(.*?).wiff", x = MM_dates$file, replacement = "\\1")
#MM_dates$ctime = strptime(as.character(MM_dates$ctime), format = "%Y-%m-%d %H:%M:%S")

#data.raw.f %>% group_by(ID, uniProt_ID, StrippedSequence, Sample, batch) %>% mutate(n = n()) %>% ungroup %>% arrange(-n,uniProt_ID, StrippedSequence, Sample, batch)
#data.raw.f.f = data.raw.f %>% group_by(ID, uniProt_ID, StrippedSequence) %>% mutate(has.na = ifelse(any(is.na(NormArea)),1,0)) %>% filter(has.na != 1) 
data.raw.f.f = data.raw.f
data.matrix = dcast(data.raw.f.f, formula=uniProt_ID+ID ~ Sample+batch, value.var = "NormArea")


data.raw.f.f$genotype[grepl(pattern = "mix", ignore.case = T, x = data.raw.f.f$Sample)] <- "QC"
MM_sample_code$Gene <- as.character(MM_sample_code$Gene)
MM_sample_code$Gene[grepl(pattern = "mix", ignore.case = T, x = MM_sample_code$Sample)] <- "QC"
MM_sample_code$Gene[MM_sample_code$SampleID == "C09b"] <- "LST7"
data.raw.f.f$genotype <- as.character(MM_sample_code$Gene[match(data.raw.f.f$Sample, MM_sample_code$Sample)])

write.table(data.raw.f.f, file = "./R/objects/data_raw.tsv", quote = F, sep = "\t", row.names = F)


data.matrix.f = tbl_df(na.omit(data.matrix))

#View(data.matrix[,c(1,2,  grep("mix", ignore.case = T, colnames(data.matrix)))])

data.matrix.f.long = data.matrix.f %>% 
  gather(Sample, NormArea, -uniProt_ID, -ID) 

data.matrix.f.long2 <- data.matrix.f.long %>% 
  separate(Sample, c("Sample", "batch"), -3) %>%
  mutate(batch = gsub("_", "", batch)) %>%
  arrange(uniProt_ID, ID)


data.matrix.f.long2.summary = data.matrix.f.long2 %>% 
  group_by(Sample, batch) %>% 
  summarize(total_signal = sum(NormArea)) %>%
  mutate(is.QC = grepl(pattern = "mix", ignore.case = T, Sample))


data.matrix.f.long2.summary = merge(data.matrix.f.long2.summary, MM_sample_code %>% filter(Batch !=1 ), by.x = c("Sample", "batch") , by.y = c("Sample", "Batch"), all = T)

#data.matrix.f.long2.summary$ctime = MM_dates$ctime[match(data.matrix.f.long2.summary$Sample, MM_dates$Sample)]
data.matrix.f.long2.summary$ctime = as.POSIXct(data.matrix.f.long2.summary$ctime)


toPlot = data.matrix.f.long2.summary

plot.before = ggplot(toPlot, aes(x=ctime, y=total_signal, shape=factor(batch), color = is.QC)) + 
  geom_point() + 
  scale_x_datetime(date_breaks = "1 day", date_labels = "%m-%d") 
plot.before
file_name = paste(fun_name,"batch_effects", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=plot.before, height=8.27, width = 2*8.27)



data.matrix.f.long3 = data.matrix.f.long2 %>%
  unite(id, uniProt_ID, ID, sep = ":") %>%
  unite(Sample_batch, Sample, batch, sep = "|")

data.matrix.f.wide = data.matrix.f.long3 %>%
  spread(Sample_batch, NormArea)

data.matrix = as.matrix(data.matrix.f.wide[,-1])
rownames(data.matrix) = data.matrix.f.wide$id

exp_metadata = as.data.frame(do.call(rbind, stringr::str_match_all(string=unique(data.matrix.f.long3$Sample_batch), pattern = "(\\w+)_SWATH_(\\d+)\\|(\\d+)")))
names(exp_metadata) = c("Sample_batch", "treatment", "replicate", "batch")
exp_metadata
pheno = exp_metadata[match(colnames(data.matrix), exp_metadata$Sample_batch),]

pheno$treatment2 = paste(pheno$treatment, pheno$batch, sep=".")
pheno$treatment2[grepl(pattern = "mix", x=pheno$treatment2)] = "mix"
      
mod = model.matrix(~as.factor(treatment2), data=pheno)

data.matrix.combat = ComBat(data.matrix, batch=pheno$batch, mod=mod, par.prior=T)
data.matrix.combat.log = ComBat(log(data.matrix), batch=pheno$batch, mod=mod, par.prior=T)


#pca stuff
before = data.matrix
after = exp(data.matrix.combat.log)

pca = prcomp(t(before), scale.=T)
x.n = 1
y.n = 2
x_var = round(pca$sdev[x.n]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[y.n]^2/sum(pca$sdev^2)*100,2)
annot = data.frame(x_var, y_var, type="before")

scores = as.data.frame(pca$x[,1:5])
scores$type = "before"
scores$sample.id = rownames(scores)

pca = prcomp(t(after), scale.=T)
x.n = 1
y.n = 2
x_var = round(pca$sdev[x.n]^2/sum(pca$sdev^2)*100,2)
y_var = round(pca$sdev[y.n]^2/sum(pca$sdev^2)*100,2)
annot = rbind(annot,data.frame(x_var, y_var, type="after"))

scores = rbind(scores, data.frame(sample.id = rownames(scores), pca$x[,1:5], type = "after"))
#scores$batch_kmeans = factor(pheno$batch[match(scores$sample.id, pheno$Sample_batch)])
scores$batch = factor(pheno$batch[match(scores$sample.id, pheno$Sample_batch)])
scores.mix = scores[grepl(pattern="mix", ignore.case=T, x=rownames(scores)),]
scores.mix$type = factor(scores.mix$type, levels=c("before", "after"))
annot$text = paste(annot$x_var, annot$y_var)

scores$type = factor(scores$type, levels = c("before", "after"))
library(cowplot)

plot.pca = ggplot(scores, aes(x=PC1, y=PC2)) + 
  geom_point(size=3, aes(col=batch) )+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(data=scores.mix, aes(x=PC1, y=PC2),size=3,col="black", shape=17) +
  geom_text(data = annot, aes(x=-50, y=-50, label=text)) +
  facet_wrap(~type, scales="fixed") + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = rel(1)))

data.matrix.combat.df = cbind(data.frame(id = rownames(data.matrix.combat.log)),
                                   as.data.frame(exp(data.matrix.combat.log)))


data.matrix.combat.long = tbl_df(data.matrix.combat.df)%>%
  gather(Sample_batch, NormArea, -id) %>%
  separate(col = "id", into = c("uniProt_ID", "ID"), sep = ":") %>%
  separate(col = "Sample_batch" , into = c("Sample", "batch"), sep="\\|")


data.matrix.combat.long$is.Corrected = 1
data.matrix.f.long2$is.Corrected = 0

final.data = rbind(data.matrix.f.long2, data.matrix.combat.long)
final.data$is.Corrected = factor(final.data$is.Corrected)

file_name = paste("final.data", fun_name, "RData", sep=".")
file_path = paste(output_dir, file_name, sep="/")
save(final.data,file=file_path)



final.data %>% filter(is.Corrected  == 1) %>% 
  group_by(uniProt_ID, Sample) %>%
  summarise(signal = exp(mean(log(NormArea))))

toPlot = final.data
plot.before = ggplot(toPlot, aes(x=ctime, y=total_signal, shape=factor(batch), color = is.QC)) + 
  geom_point() + 
  scale_x_datetime(date_breaks = "1 day", date_labels = "%m-%d") 

data.final.summary = final.data %>% group_by(Sample, batch, is.Corrected) %>% summarize(total_signal = sum(NormArea))
data.final.summary$is.QC = grepl(pattern = "mix", ignore.case = T, x = data.final.summary$Sample)


data.final.summary = merge(data.final.summary, MM_sample_code, by.x = c("Sample", "batch") , by.y = c("Sample", "Batch"))
#data.final.summary$ctime = MM_dates$ctime[match(data.final.summary$Sample, MM_dates$Sample)]
data.final.summary$ctime = as.POSIXct(data.final.summary$ctime)


file_name = paste("data.final.summary", fun_name, "RData", sep=".")
file_path = paste(output_dir, file_name, sep="/")
save(data.final.summary,file=file_path)


#MM_sample_code %>% filter(SampleID == "mix") %>% group_by(Batch) %>% summarise(n= n())
#MM_sample_code %>% arrange(ctime) %>% filter(Batch != 1)

toPlot = data.final.summary




p = ggplot(toPlot, aes(x=ctime, y=total_signal, shape=factor(batch), color = is.QC)) + 
  geom_point() + 
  scale_x_datetime(date_breaks = "1 day", date_labels = "%m-%d") +
  facet_grid(~is.Corrected) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        aspect.ratio = 3/8)

p.final = plot_grid(p, plot.pca, labels = c("A", "B"), nrow = 2, align = "v")

file_name = paste(fun_name,"batch_effects", "pdf", sep=".")
file_path = paste(figures_dir, file_name, sep="/")
ggsave(filename=file_path, plot=p.final, height=8.27, width = 2*8.27)



