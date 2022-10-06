library(readr)
library(ggseqlogo)
library(dplyr)
library(Biostrings)
library(ggplot2)

# Import file
file1 <-list.files(path = ".", pattern = "filtered_final.txt")
# convert to table
table.k<-read.table(file1, header = T)

#setwd("/Users/lucasmatiasferreri/OneDrive - Emory University/Collaborations/Nahara/Playpen01_Mitigate/NGS-Barcode/DataAnalysis/Serial_dilutions/NV-024/NV-024_01")

#table.k<- read.table("filtered_final.txt", sep = "\t", header = T)

##Convert sequence to DNAString
table.k$nt.DNAString<- DNAStringSet(table.k$nt.DNAString)

table.k[c("Codons")] <- NA
# Parse codons
table.k$Codons <- paste(
  subseq(table.k$nt.DNAString, 41, 43, NA),
  subseq(table.k$nt.DNAString, 47, 49, NA),
  subseq(table.k$nt.DNAString, 53, 55, NA),
  subseq(table.k$nt.DNAString, 59, 61, NA),
  subseq(table.k$nt.DNAString, 62, 64, NA),
  subseq(table.k$nt.DNAString, 68, 70, NA),
  subseq(table.k$nt.DNAString, 71, 73, NA),
  subseq(table.k$nt.DNAString, 74, 76, NA),
  subseq(table.k$nt.DNAString, 80, 82, NA),
  subseq(table.k$nt.DNAString, 83, 85, NA),
  subseq(table.k$nt.DNAString, 98, 100, NA),
  subseq(table.k$nt.DNAString, 101, 103, NA),
  sep = ""
)

# AA column
table.k[c("AA")]<- NA
# Convert to StringSet
table.k$Codons<-DNAStringSet(table.k$Codons)
# Translate
table.k$AA <- translate(table.k$Codons, if.fuzzy.codon = "X", genetic.code = GENETIC_CODE)

table.k$nt.DNAString <- as.character(table.k$nt.DNAString)
table.k$Codons <- as.character(table.k$Codons)
table.k$AA <- as.character(table.k$AA)

#Calculate frequency and filter
multi.aa <- table.k %>%
  group_by(AA) %>%
  summarise(n = n()) %>% 
  unique() %>% 
  mutate(freq = n / sum(n)) 

# Add Sample and Index
y <- table.k %>% 
  select(Sample, Index, AA) %>% 
  unique()

multi.aa <-
multi.aa %>% 
  left_join(y, by = "AA")

write.table(multi.aa, 
            file = "multi.aa.txt", 
            quote = F, sep="\t",
            col.names = F, 
            row.names = F)

# create logos plot
logo <- 
  ggseqlogo(table.k$AA, method = 'prob', facet = "wrap", scales = "free_x", ncol = NULL, nrow = NULL) 
ggsave("logos.plot.pdf", dpi = 350, height = 3, width = 6)
