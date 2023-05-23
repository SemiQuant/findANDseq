#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
require(openPrimeR, quietly = T)
require(stringr, quietly = T)

seq.df <- read_templates(args[1])
# seq.df <- seq.df[seq.df$ID != '>MTB000019|rrs|rRNA|1471846-1473382|+|Ribosomal RNA 16S', ]
# seq.df <- seq.df[seq.df$ID != '>MTB000020|rrl|rRNA|1473658-1476795|+|Ribosomal RNA 23S', ]
# seq.df <- seq.df[grepl("CDS", seq.df$ID),]
seq.df <- seq.df[!grepl("rRNA|Ribosomal RNA", seq.df$ID, ignore.case = T), ]
# seq.df$ID <- gsub("\\|.*", "", seq.df$ID)
seq.df$ID <- str_replace_all(seq.df$ID, "[[:punct:]]", " ")
seq.df$ID <- strtrim(seq.df$ID, 50)# to stop filename too long
seq.df$Allowed_Start_rev <- seq.df$Allowed_Start_rev_ali <- seq.df$Allowed_Start_rev_initial <- seq.df$Allowed_Start_rev_initial_ali <- ifelse(seq.df$Sequence_Length >= 500, 200,  #ifelse(seq.df$Sequence_Length < 200, 50, seq.df$Sequence_Length-150)
                                                                                                                                            ifelse(seq.df$Sequence_Length < 100, 50, 
                                                                                                                                                ifelse(seq.df$Sequence_Length < 200, 100, 
                                                                                                                                                    50
                                                                                                                                                    )))

seq.df$Allowed_End_rev <- seq.df$Allowed_End_rev_ali <- seq.df$Allowed_End_rev_initial <- seq.df$Allowed_End_rev_initial_ali <-  ifelse(seq.df$Sequence_Length < 500, seq.df$Sequence_Length, 500) #seq.df$Sequence_Length-80

if (nrow(seq.df[seq.df$Sequence_Length < 80, ] > 0)){
  print("Removing the following as <80nts: ")
  print(seq.df[seq.df$Sequence_Length < 80, ]$ID)
  seq.df <- seq.df[seq.df$Sequence_Length >= 80, ]
}
# seq.df$region <- paste0(seq.df$Allowed_Start_rev, ',', seq.df$Allowed_End_rev-seq.df$Allowed_Start_rev)
seq.df_out <- data.frame(ID = gsub('>', '', seq.df$ID),
                         Sequence = seq.df$Sequence,
                         Start = seq.df$Allowed_Start_rev,
                         Length = seq.df$Allowed_End_rev-seq.df$Allowed_Start_rev
                         # region = seq.df$region
                         )
write.table(seq.df_out, args[2], sep = "\t", quote = F, row.names = F, col.names = T)