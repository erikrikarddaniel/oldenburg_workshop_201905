##########Atacama Qiime2 example 
#after running the qiime2biom2tsv.md file, 
#you can import the data like this as well.


###load libraries

library(readr)
library(dplyr)
library(tidyr)

asvs <- read_tsv(
  'feature-table.tsv', 
  col_types = cols(.default = col_double(), `#OTU ID` = col_character()), 
  skip = 1
) %>% 
  rename(asv_id = `#OTU ID`) %>% 
  gather(sample, count, 2:ncol(.)) %>% 
  filter(count > 0) %>% 
  mutate(count = as.integer(count))

#unzip the taxonomy file taxonomy.qza
#get the taxomony
taxonomy <- read_tsv( 'taxonomy.tsv')  %>% 
    rename(asv_id =`Feature ID`)



#get metadata
metadata <- read_tsv('metadata.yml')


