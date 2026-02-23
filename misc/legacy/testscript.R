library(tidyverse)
library(devtools)
load_all()

df <- read.table(file="../../inst/extdata/WT_6mA_Mg.txt", header=T)
meta <- read_tsv('../../inst/extdata/all_site_annotations.txt')

df_annotated <- annotateMethylSites(df, meta, location='Position')
