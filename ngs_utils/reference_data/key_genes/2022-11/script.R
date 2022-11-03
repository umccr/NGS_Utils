require(tidyverse)
require(here)
require(glue)

wd <- here("ngs_utils/reference_data/key_genes/2022-11")
cancermine <- file.path(wd, "sources/cancermine_collated.tsv") |>
  read_tsv(col_types = cols(.default = "c", citation_count = "i"))
cosmic <- bind_rows(
  file.path(wd, "sources/Census_all_Tier1.tsv") |> read_tsv(),
  file.path(wd, "sources/Census_all_Tier2.tsv") |> read_tsv())
cpsr <- file.path(wd, "sources/cpsr_superpanel_2022_01.xlsx") |>
  readxl::read_xlsx()

ncg <- file.path(wd, "sources/NCG_cancerdrivers_annotation_supporting_evidence.tsv") |>
  read_tsv(col_types = cols(.default = "c"))

# oncokb differences
# 2022-10-28, 1,068 genes
oncokb <- file.path(wd, "sources/cancerGeneList.tsv") |>
  read_tsv(col_types = cols(.default = "c"))
oncokb_old <- file.path(wd, "../sources/oncoKB_cancerGeneList.txt") |>
  read_tsv(col_types = cols(.default = "c"))

table(oncokb$`Hugo Symbol` %in% oncokb_old$`Hugo Symbol`)
table(oncokb_old$`Hugo Symbol` %in% oncokb$`Hugo Symbol`)
oncokb[!oncokb$`Hugo Symbol` %in% oncokb_old$`Hugo Symbol`, ] |> arrange(`Hugo Symbol`)
oncokb_old[!oncokb_old$`Hugo Symbol` %in% oncokb$`Hugo Symbol`, ] |> arrange(`Hugo Symbol`)

# TEMPUS differences
tempus <- file.path(wd, "sources/tempus.txt") |>
  read_tsv(col_types = "c", col_names = "symbol") |>
  mutate(symbol = sub("\\*+", "", symbol), # remove stars
         symbol = sub("\\(", "", symbol), # remove parentheses
         symbol = sub("\\)", "", symbol))
tempus_old <- file.path(wd, "../sources/TEMPUS.genes") |>
  read_tsv(col_types = "c", col_names = "symbol")

table(tempus$symbol %in% tempus_old$symbol)
table(tempus_old$symbol %in% tempus$symbol)

# old ones not in new
tempus_old[!tempus_old$symbol %in% tempus$symbol, ]
# new ones not in old
tempus[!tempus$symbol %in% tempus_old$symbol, ]
