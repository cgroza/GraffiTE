library(repeatR)
library(dplyr)
library(readr)
library(tidyverse)
library(RColorBrewer)

te_df <- read_rm("Zea_bz_lr.fasta.out", include_secondary = T, tibble = T) %>%
  mutate(aID = row_number())

node_sizes <- read_delim("node_sizes.csv", delim = " ", col_names = c("node", "size"), col_types = "ci")

path_df <- read_tsv("Zea_bz_lr.gaf", col_names = F) %>%
  select(X1, X6) %>%
  dplyr::rename(qname = X1, path = X6)

df <- left_join(path_df, te_df, by = "qname") %>%
  mutate(path = str_sub(path, 2)) %>%
  mutate(node = path) %>%
  group_by(aID) %>%
  separate_rows(node, sep = "[><]")

df <- left_join(df, node_sizes) %>%
  group_by(path, aID) %>%
  mutate(
    pend = cumsum(size),
    pstart = pend - size + 1,
    ) %>%
  mutate(
    tname = case_when(
      pstart >= qstart & pend <= qend ~ tname,
      T ~ NA
    ),
    tclass = case_when(
      pstart >= qstart & pend <= qend ~ tclass,
      T ~ NA
    ),
    ) %>%
  ungroup()

te_families <- unique(na.omit(df$tclass))
te_names <- unique(na.omit(df$tname))
te_names <- te_names[!str_detect(te_names, "chr[49]_")]

color_families_df <- tibble(
  tclass = te_families,
  Color = colorRampPalette(brewer.pal(12, "Paired"))(length(te_families))
)

color_names_df <- tibble(tname = te_names,
                         Color = brewer.pal(length(te_names), "Paired"))

left_join(df, color_names_df, by = "tname") %>%
  select(node, tclass, tname, Color) %>%
  filter(!is.na(tclass)) %>%
  filter(!is.na(Color)) %>%
  dplyr::rename(Name = node) %>%
  write_csv("bz_bandage_name_test.csv")

left_join(df, color_families_df, by = "tclass") %>%
  select(node, tclass, tname, Color) %>%
  filter(!is.na(tclass)) %>%
  filter(!is.na(Color)) %>%
  dplyr::rename(Name = node) %>%
  write_csv("bz_bandage_family.csv")
