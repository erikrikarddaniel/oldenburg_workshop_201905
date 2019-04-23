### Run dada2 in R

#install package from bioconductor if required
#Note, if you are using a Mac (like me), you might have to update R
#to a version >3.5, otherwise dada2 might give errors while installing.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")

library(dada2)
#library(vegan)
library(ggplot2)
library(tidyverse)
packageVersion("dada2")

#citation("dada2") #for credits
#The below script is modified from the RaukR course
#https://nbisweden.github.io/RaukR-2018/metagenomics_John/lab/DADA2_workshop.html

#set working direktory
setwd("/Users/cbunse/Desktop/CoDa/gitpage/testscripts/MiSeq_Atacama0.01/qiime2-atacama-tutorial/")
#your data need to be demultiplexed, 
#                     primers/adapters etc removed, 
#                     if paired seqs then in matched order

###################
#since this is not the case for the Atacama-desert-data example (with 1% of the data), 
#we have to demultiplex them first and we will do that in the terminal/qiime2 way,
# see https://docs.qiime2.org/2019.1/tutorials/atacama-soils/#

  #load qiime $source activate qiime2-2019.1
  #install the atacama dataset and metafile, rename and store in right folder
      #(see qiime2 tutorial)
  #create EMPPairedEndSequences artefact
  #demultiplex
  #summarize qiime demux file

  #export demux data to fastqfiles

###################

# Now in R
#load and check out the metafile

dataframe <- read.csv("sample-metadata.txt", sep="\t", head=TRUE)
head(dataframe) # you will realize that the second row does not contain data, so let's remove it
df<- dplyr::slice(dataframe, 3:n())
head(df)

#load the demultiplexed data!!
demux<- read.csv("demux.qza", sep="\t", head=FALSE)


#store the names of the sample names and extract them
path <- paste(getwd(), "data", sep="/")
fR1 <- sort(list.files(path, pattern="L1_R1.fastq.gz", full.names = TRUE))
fR2 <- sort(list.files(path, pattern="L1_R2.fastq.gz", full.names = TRUE))
sample.names <- regmatches(basename(fR1),regexpr("[0-9]+",basename(fR1)))


#trimm and filter
plotQualityProfile(c(fR1[1],fR2[1]))

# specify the output directory and file names for the filtered reads:
filt_path <- file.path(path, "filtered")
filtfR1 <- file.path(filt_path, paste0(sample.names, "_R1_filt.fastq.gz"))
filtfR2 <- file.path(filt_path, paste0(sample.names, "_R2_filt.fastq.gz"))

filt_trim <- filterAndTrim(fR1, filtfR1, fR2, filtfR2, truncLen=c(200,150), trimLeft = 20,
                           maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                           compress=TRUE, multithread=TRUE)

#dereplication
derepR1 <- derepFastq(filtfR1, verbose=TRUE)
names(derepR1) <- sample.names
derepR2 <- derepFastq(filtfR2, verbose=TRUE)
names(derepR2) <- sample.names

#estimate error model
error_ratesR1 <- learnErrors(derepR1, multithread=TRUE)
error_ratesR2 <- learnErrors(derepR2, multithread=TRUE)

plotErrors(error_ratesR1)
plotErrors(error_ratesR2)

#get ASVs
dadaR1 <- dada(derepR1, err=error_ratesR1, multithread=TRUE)
dadaR2 <- dada(derepR2, err=error_ratesR2, multithread=TRUE)

#merge
mergers <- mergePairs(dadaR1, derepR1, dadaR2, derepR2, verbose=TRUE)

#seq table
seqtab <- makeSequenceTable(mergers)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE)

#check chimeras & proportions
print(1-dim(seqtab.nochim)[2] / dim(seqtab)[2])
print((1-sum(seqtab.nochim)/sum(seqtab))*100)

#check sequences throughout
getN <- function(x) sum(getUniques(x))
track <- cbind(filt_trim, sapply(dadaR1, getN), sapply(dadaR2, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
datatable(track, options = list(pageLength=24))

#taxonomy
system("curl -o data/rdp_species_assignment_16.fa.gz https://zenodo.org/record/801828/files/rdp_species_assignment_16.fa.gz")
system("curl -o data/rdp_train_set_16.fa.gz https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz")

taxa <- assignTaxonomy(seqtab.nochim, "data/rdp_train_set_16.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "data/rdp_species_assignment_16.fa.gz")






