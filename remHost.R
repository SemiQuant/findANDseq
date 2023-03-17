#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
require(tidyverse)
require(Biostrings)

host <- read_tsv(args[1], col_names = F)[c(1)]
file_in <- args[2]
gene_names <- read_tsv(args[3])$ID

dat <- read_tsv(file_in, col_names = F, col_select = c(1,3,4,5), col_types = list("c", "c", "f", "i", "f"))
colnames(dat) <- c("ID", "Align", "Start", "Seq")

# dat$Align <- as.factor(gsub("\\|.*", "", dat$Align))

host <- host %>% 
  separate(X1, c("host_genes", "ID"), sep = " ") %>% 
  mutate(host_genes = as.numeric(host_genes))

dat <- dat %>%
  filter(Start >= 50) %>% 
  filter(!ID %in% host[host$host_genes > 1,]$ID) %>% 
  filter(as.character(Align) %in% gene_names) 

write_tsv(dat, file = args[4])
