# Export biom files QIIME2's artefacts and convert to tab separated (tsv) files

## Export a biom file from a QIIME2 artefact

TODO

## Install `biom-format` with Conda

```{bash}
conda create -n biom-format
conda activate biom-format
conda install -c conda-forge biom-format
conda install -c conda-forge h5py
```

## Converting to tsv and importing into R with Tidyverse

After installing the biom-format package, you can convert biom files to tsv.

```{bash}
conda activate biom-format
biom convert -i feature-table.biom --to-tsv -o feature-table.tsv
```

Read the tab separated file into R using the `readr` package's `read_tsv` function,
and do some tidying up with functions from `dplyr` and `tidyr`.

```{R}
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
```
