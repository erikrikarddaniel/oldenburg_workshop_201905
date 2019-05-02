
####### CoDa approach for Amplicon data #######

library(CoDaSeq) # as always, check the vignette if questions arise
library(ALDEx2)
library(zCompositions)
library(tidyverse)

#data(ak_op)
#data(hmpgenera)

#read-in output from "qiime2biom2tsv-script"
asvs <- read_tsv(
  'feature-table.tsv', 
  col_types = cols(.default = col_double(), `#OTU ID` = col_character()), 
  skip = 1
) %>% 
  rename(asv_id = `#OTU ID`) %>% 
  gather(sample, count, 2:ncol(.)) %>% 
  filter(count > 0) %>% 
  mutate(count = as.integer(count))

#convert to wide format so that it is conform to downstream analyses
wasvs<- asvs %>%
  tidyr::spread(sample, count, fill=0)%>%
  tibble::column_to_rownames("asv_id")

## optional: subset and filter the dataset to a certain amount of reads, proportion, ...
f <- codaSeq.filter(wasvs, min.reads=10, min.prop=0.001, max.prop=1,
                    min.occurrence=0.25,  samples.by.row=FALSE)

## replace zeros (as we will deal with geomeans later)
# also, transpose the datafile
f.n0 <- zCompositions::cmultRepl(t(f), label=0, method="CZM")

# perform clr transformation
f.clr <- codaSeq.clr(f.n0, samples.by.row=TRUE)


# perform the singular value decomposition on clr transformed values
f.pcx <- prcomp(f.clr)
biplot(f.pcx, var.axes=FALSE, cex=c(0.5,0.6), scale=0)


### Aitchicon Distance
m<- as.matrix(compositions::dist(f.clr))
m[upper.tri(m, diag = FALSE)] <- NA #remove redundant comparisons

#example tidy data for ggplot
mw <- reshape::melt(m)  %>%  # bring data into right format
      na.omit()   # exclude NAs
names(mw) <- c('row', 'column', 'value')  # tidy data
mw$row <- as.character(mw$row)  # tidy data
mw$column <- as.character(mw$column)  # tidy data

# example plot
ggplot(mw, aes(row, column)) + 
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient(high = "white",low = "steelblue")+ 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  ggtitle("Aitchison Distance, Atacama example")+
  theme(axis.text.x = element_text(size=6, angle=45, hjust = 1))+
  theme(axis.text.y = element_text(size=6))              


# Ordination (OCA of variance)
#singular value decomposition
x<- prcomp(f.clr)
biplot(x)

##################################################################

####Stripchart example if you have different treatments, run the above commands for the ak_op dataset, then
# perform the singular value decomposition on clr transformed values
conds <- c(rep("A", 15), rep("O", 15))
f.x <- aldex.clr(f, conds)
f.e <- aldex.effect(f.x, conds)
f.t <- aldex.ttest(f.x, conds)
f.all <- data.frame(f.e,f.t)

codaSeq.stripchart(aldex.out=f.all, group.table=hmpgenera, group.label="genus",
                   sig.cutoff=0.05, sig.method="we.eBH", x.axis="effect",
                   mar=c(4,8,4,0.5))

######differential abundance

#ALDEx2 example (see vignette for more infos), note that ALDEx2 calculates the clr values itself
data(selex)
#subset for efficiency of example
selex <- selex[1201:1600,]
#define experimental conditions / groups of samples
conds <- c(rep("NS", 7), rep("S", 7))

x <- aldex(selex, conds, mc.samples=16, test="t", effect=TRUE,
          include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)
aldex.plot(x, type="MA", test="welch")
aldex.plot(x, type="MW", test="welch")

#clr transformation in ALDEx2
x <- aldex.clr(selex, conds, mc.samples=16, denom="iqlr", verbose=TRUE)



