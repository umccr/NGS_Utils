require(biomaRt)
require(tidyverse)

umccr_genes <- "umccr_cancer_genes.latest.genes" |>
  readr::read_tsv(col_names = "gene", col_types = "c") |>
  dplyr::pull(gene)
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 106)
attributes <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position")
filters <- "hgnc_symbol"
values <- list("hgnc_symbol" = umccr_genes)
d <- getBM(attributes=attributes, filters=filters, values=values, mart=mart) |>
  as_tibble()

umccr_unknown <- umccr_genes[!umccr_genes %in% d$hgnc_symbol]

# Let's manually grab the HGNC symbols..
# First the ones found from the Ensembl website
unknown1 <- tibble::tribble(
  ~gene,  ~alias,
  "C2orf40", "ECRG4",
  "C2orf44", "WDCP",
  "CARS", "CARS1",
  "CASC5", "KNL1",
  "CYR61", "CCN1",
  "FAM175A", "ABRAXAS1",
  "FAM46C", "TENT5C",
  "FAM58A", "CCNQ",
  "FGFR1OP", "CEP43",
  "H3F3A", "H3-3A",
  "H3F3B", "H3-3B",
  "H3F3C", "H3-5",
  "HIST1H1B", "H1-5",
  "HIST1H1C", "H1-2",
  "HIST1H1D", "H1-3",
  "HIST1H1E", "H1-4",
  "HIST1H2AC", "H2AC6",
  "HIST1H2AG", "H2AC11",
  "HIST1H2AL", "H2AC16",
  "HIST1H2AM", "H2AC17",
  "HIST1H2BC", "H2BC4",
  "HIST1H2BD", "H2BC5",
  "HIST1H2BG", "H2BC8",
  "HIST1H2BJ", "H2BC11",
  "HIST1H2BK", "H2BC12",
  "HIST1H2BO", "H2BC17",
  "HIST1H3A", "H3C1",
  "HIST1H3B", "H3C2",
  "HIST1H3C", "H3C3",
  "HIST1H3D", "H3C4",
  "HIST1H3E", "H3C6",
  "HIST1H3F", "H3C7",
  "HIST1H3G", "H3C8",
  "HIST1H3H", "H3C10",
  "HIST1H3I", "H3C11",
  "HIST1H3J", "H3C12",
  "HIST1H4I", "H4C9",
  "HIST2H3A", "H3C15",
  "HIST2H3C", "H3C14",
  "HIST2H3D", "H3C13",
  "HIST3H3", "H3-4",
  "KIAA1598", "SHTN1",
  "LHFP", "LHFPL6",
  "MKL1", "MRTFA",
  "MLLT4", "AFDN",
  "MRE11A", "MRE11",
  "PAK7", "PAK5",
  "PARK2", "PRKN",
  "PDL1", "CD274",
  "PDL2", "PDCD1LG2",
  "RFWD2", "COP1",
  "SEPT5", "SEPTIN5",
  "SEPT6", "SEPTIN6",
  "SEPT9", "SEPTIN9",
  "TAZ", "TAFAZZIN",
  "TCEB1", "ELOC",
  "TNFRSF6", "FAS",
  "WHSC1", "NSD2",
  "WHSC1L1", "NSD3",
  "WISP1", "CCN4",
  "WISP3", "CCN6",
  "ZNF198", "ZMYM2",
  "ZNF278", "PATZ1"
)

# Next the ones found via GeneCards
unknown2 <- tibble::tribble(
  ~gene,  ~chr, ~start, ~end,
  "CCAT1", "8", 127207381, 127219268,
  "IGH", "14", 105586437, 106879844,
  "IGK", "2", 88857361, 90235368,
  "IGL", "22", 22026076, 22922913,
  "TRA", "14", 21621904, 22552132,
  "TRB", "7", 142299011, 142813287,
  "TRD", "14", 22422546, 22466577
)

values2 <- list("hgnc_symbol" = unknown1$alias)
d2 <- getBM(attributes=attributes, filters=filters, values=values2, mart=mart) |>
  as_tibble()
d3 <- bind_rows(d, d2)

# 1,232 genes
x1 <- d3 |>
  filter(!grepl("_", chromosome_name)) |>
  distinct(hgnc_symbol) |>
  pull(hgnc_symbol)

# which umccr genes don't we have?
x2 <- umccr_genes[!umccr_genes %in% x1]
# those should be in the gene column of unknown1
table(x2 %in% unknown1$gene)
# and the rest should be in unknown2
x2[!x2 %in% unknown1$gene]
x2[!x2 %in% unknown1$gene] %in% unknown2$gene

# Now gather all
res <- d3 |>
  filter(!grepl("_", chromosome_name)) |>
  select(gene = hgnc_symbol, chr = chromosome_name, start = start_position, end = end_position) |>
  bind_rows(unknown2) |>
  # mutate(chr = paste0("chr", chr)) |>
  distinct(gene, .keep_all = TRUE) |> # RMRP is duplicated, just keep the first one
  arrange(gene)

# add +1000bp buffer just to be safe
pad <- 1000
res <- res |>
  mutate(start = start - pad, end = end + pad) |>
  select(chr, start, end) |>
  arrange()
