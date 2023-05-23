#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
require(tidyverse, quietly = T)
require(Biostrings, quietly = T)

host <- read_tsv(args[1], col_names = F)[c(1)]
file_in <- args[2]
gene_names <- read_tsv(args[3])$ID

dat <- read_tsv(file_in, col_names = F, col_select = c(1,3,4,5), col_types = list("c", "c", "f", "i", "f"))
colnames(dat) <- c("ID", "Align", "Start", "Seq")

host <- host %>% 
  separate(X1, c("host_genes", "ID"), sep = " ") %>% 
  mutate(host_genes = as.numeric(host_genes))

dat <- dat %>%
  filter(ID %in% host$ID)
  
# change this to check if removing all of one gene primerset
dat_red <- dat %>%
  filter(Start >= 50) %>%
  filter(!ID %in% host[host$host_genes > 1,]$ID) %>% 
  filter(as.character(Align) %in% gene_names) 

# not so goods
host <- host %>% 
  arrange(desc(host_genes)) %>% 
  filter(host_genes >= min(host_genes) - 2 & host_genes <= min(host_genes) + 2)
  # filter(row_number() <= 50)

dat_notSoGood <- dat %>% 
  filter(gsub("_.*", "", ID) %in% gene_names[!gene_names %in% gsub("_.*", "", dat_red$ID)]) %>% 
  arrange(desc(nchar(as.character(Seq))), desc(ifelse(Start >= 300 & Start <= 500, Start * 1000, Start))) %>% 
  group_by(Align) %>% 
  filter(row_number() <= 10)

dat <- rbind(dat_red,
             dat_notSoGood)

write_tsv(dat, file = args[4])
