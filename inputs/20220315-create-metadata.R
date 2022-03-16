library(readr)
library(dplyr)

sra_annotations <- read_tsv("https://raw.githubusercontent.com/greenelab/core-accessory-interactome/master/data/metadata/SRA_annotations.tsv",
                            skip = 1, col_names = c("experiment", "strain_type")) %>%
  filter(strain_type %in% c("PAO1", "PA14"))

write_tsv(sra_annotations, "inputs/metadata.tsv")
