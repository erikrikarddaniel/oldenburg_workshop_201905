##########Atacama Qiime2 example for course

###load libraries

library(tidyverse)
library(biomformat) ## I used this one instead of biom, but still had some issues installing it
library(phyloseq)
library(breakaway)

##import biom dataset
dat <- read_biom("feature-table.biom")

#Coerce dat into a matrix and then into a data frame:
asv_table <- as.data.frame(as.matrix(biom_data(dat)))
#get the taxomony
taxonomy <- observation_metadata(dat)


#get metadata
metadata <- sample_metadata(dat)
## not available in atacama example?


############alpha diversity tests (from breakaway example)
asv_m<- asv_table %>% 
  #dplyr::select(2:6) %>%
  remove_rownames()

head(asv_m)

breakaway(asv_m)
breakaway(apples)

