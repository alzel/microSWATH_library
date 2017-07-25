library("RColorBrewer") ## Color palettes
library("ggplot2") ## Convenient and nice plotting
library("reshape2") ## Flexibly res
library("mzR")
library("dplyr")
#library(MSnbase)
library(xcms)

rm(list = ls())

library.raw <- read.delim("./results/2015-12-05/SpecLib_cons_openswath_16.csv")
library.filtered = library.raw[grep(pattern = "Biognosys",ignore.case = T, x = library.raw$ProteinName),]

file = "results/2015-12-05/data/microSWATH_library_01.mzXML"
xr <- xcmsRaw(file, includeMSn=TRUE)
xcmsSet()

rtrange.m= matrix(rep(c(100,1000),nrow(library.filtered)),nrow = nrow(library.filtered), byrow = T)
mzrange.m = cbind(library.filtered$ProductMz - 0.055, + library.filtered$ProductMz + 0.055)
#mzrange.m = cbind(487.2567 - 0.055, + 487.2567 + 0.055)
rtrange.m= matrix(rep(c(3.336,1800.18),nrow(mzrange.m)),nrow = nrow(mzrange.m), byrow = T)

test = getEIC(xraw, mzrange=mzrange.m, step = 1, rtrange = rtrange.m) 
?getEIC
my_list = test@eic$xcmsRaw
xic.matrix = do.call(rbind,my_list)

dots_number=unique(unlist(lapply(my_list, nrow)))
stopifnot(length(dots_number) == 1)

transition_name = rep(library.filtered$transition_name, each=dots_number)
xic.df = as.data.frame(xic.matrix)
xic.df$transition_name = transition_name

names(xic.df) = c("rt", "intensity", "transition_name")  
names(library.raw)
annotations = dplyr::select(library.raw, transition_name, PrecursorMz, PeptideSequence)
            
xic.df.annotated = merge(xic.df, annotations, by.x = "transition_name" , by.y = "transition_name")
names(xic.df.annotated)

xic.df.summary = xic.df.annotated %>% group_by(PeptideSequence, PrecursorMz, rt) %>% summarise(intensity.sum = sum(intensity, na.rm=T)) %>% 
                                               group_by(PeptideSequence, PrecursorMz) %>% mutate(apex=rt[which.max(intensity.sum)])

ggplot(xic.df.summary, aes(x=rt, y=intensity.sum, colour=PeptideSequence)) + 
      geom_vline(aes(xintercept=apex)) + 
      geom_line()




lapply(my_list, 
       FUN = function(x) {
         print(as.character(x[[]))
       })

rep

min(test2$rt)
max(test2$rt)

?xic
str(test)
plot(e)

xic(object = ms, )

?xcmsRaw



xraw <- xcmsRaw(file, includeMSn=TRUE)
peaks <- findPeaks(xraw)

View(peaks)
extractXIC = function(object, s_peaks, mz, width = 0.055) {
  mz = 100
  width = 100
  lowMz = mz - width
  highMz = mz + width
  s_peaks = selected_peaks
  
  selected_peaks.df = lapply(s_peaks, function(x) {x[x[,1] > lowMz  & x[,1] < highMz,]})
  selected_peaks = peaks(ms,spectra_idx)
  selected_peaks.df = lapply(selected_peaks, function(x) {as.data.frame(x)})
  
  names(selected_peaks.df) = spectra_idx
  
  selected_peaks.df = selected_peaks.df[lapply(selected_peaks.df, nrow)>0]
  
  tmpnames = as.numeric(as.character(names(selected_peaks.df)))
  
  for (j in 1:length(selected_peaks.df)) {
      selected_peaks.df[[j]]$name = tmpnames[j]  
  }
  
  selected_peaks.matrix = lapply(selected_peaks.df, function(x) {as.matrix(x)})
  my_peaks = do.call(rbind, selected_peaks.matrix)
  my_peaks = as.data.frame(my_peaks)
  my_peaks$selected.ion = 100
  my_peaks$width = 100
  
  colnames(my_peaks) = c("mz", "Intensity", "acquisitionNum", "selected.ion", "width")
  return(my_peaks)
}

my_peaks$rt = hd$retentionTime[match(my_peaks$acquisitionNum, hd$acquisitionNum)]
my_peaks$mz2 = round(my_peaks$mz)
my_peaks$rt2 = round(my_peaks$rt)

my_peaks.summary = my_peaks %>% group_by(mz2, rt2) %>% summarise (sumI = sum(Intensity))

my_peaks.summary.wide = dcast(my_peaks.summary, formula = mz2~rt2, value.var = "sumI" )

my_peaks.summary.wide[is.na(my_peaks.summary.wide)] = 1

library(plotly)
View(volcano )

toPlot = as.matrix(my_peaks.summary.wide)[,-1]
pdf("test.pdf", width = 11.69, height = 8.27)

x <- list(
  title = "m/z",
)
y <- list(
  title = "Time",

)

plot_ly(z= log(toPlot), type = "surface") %>% layout(xaxis = x, yaxis = y)
dev.off()
?plot_ly
View(toPlot)
test = extractXIC(object = ms, s_peaks = selected_peaks.df, mz = Z$ProductMz,  width = 0.055 )
View(test)
library(plyr)
selected_peaks
head(library.filtered)
my.tmp = ddply(library.filtered, .variables = .(transition_name, ProductMz), 
            .fun = function(x, sp=selected_peaks.df) {
              Z <<-x
              spX <<- sp
              z = extractXIC(object = ms, s_peaks = sp, mz=x$ProductMZ, width=0.055)
              return(z)
            })
asnames(spX)
Z$ProductMz
Z
View(test)
plot(test$name, test$`sum(V2)`)
library(xcms)
library(ggplot2)
file = "results/2015-12-05/data/microSWATH_library_01.mzXML"

??xcms
xraw <- xcmsRaw(file, includeMSn=TRUE)
?xi
peaks = findPeaks(xraw)
peaks.df = as.data.frame(peaks)
View(peaks.df)
ggplot(peaks.df, aes(x=rt, y = mz)) + geom_point()
View(peaks.df)


