# Instructions for how to convert QIIME2's biom files to tab separated (tsv) files

## Install `biom-format` with Conda

```{bash}
conda create -n biom-format
conda activate biom-format
conda install -c conda-forge biom-format
conda install -c conda-forge h5py
```

## Converting to tsv and importing into R with Tidyverse

```{bash}
```

```{R}
library(readr)
library(dplyr)
bt <- read_tsv('feature-table.tsv', skip = 1) %>% 
  rename(asv_id = `#OTU ID`)
```
