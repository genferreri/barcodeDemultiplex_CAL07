# This script calculates the frequency of each barcode and the Shannon entropy. 
# It also plots the frequency distribution of the barcode species present in the sample.
# 20220921

library(readr)
library(dplyr)
library(ggplot2)
library(Biostrings)
library(vegan)
library(RColorBrewer)
library(viridis)

# Import file
file1 <-list.files(path = ".", pattern = "\\.txt")
# convert to table
table.x<-read.table(file1, header = F)

#setwd("/Users/lucasmatiasferreri/OneDrive - Emory University/Collaborations/Nahara/Playpen01_Mitigate/NGS-Barcode/DataAnalysis/Serial_dilutions/NV-024")
#table.x<-read.table("NV-024_01_final.txt", header = F)

colnames(table.x) <- c("Seq","Sample","Index")

table.x[c("nt.DNAString", "motif")]<-NA

##Convert sequence to DNAString
table.x$nt.DNAString<- DNAStringSet(table.x$Seq)
# Remove original column
table.x<-table.x[,-1]

## Parse the sites in the barcode for Pan99
#table.x$motif <- paste(
#  subseq(table.x$nt.DNAString, 51, 51, NA),
#  subseq(table.x$nt.DNAString, 54, 54, NA),
#  subseq(table.x$nt.DNAString, 58, 58, NA),
#  subseq(table.x$nt.DNAString, 72, 72, NA),
#  subseq(table.x$nt.DNAString, 75, 75, NA),
#  subseq(table.x$nt.DNAString, 81, 81, NA),
#  subseq(table.x$nt.DNAString, 87, 87, NA),
#  subseq(table.x$nt.DNAString, 90, 90, NA),
#  subseq(table.x$nt.DNAString, 91, 91, NA),
#  subseq(table.x$nt.DNAString, 93, 93, NA),
#  subseq(table.x$nt.DNAString, 96, 96, NA),
#  subseq(table.x$nt.DNAString, 99, 99, NA),
#  sep = ""
#)

## Parse the sites in the barcode for Cal07. 
#the sequence below correspond to numbering considering the 4 nt from the index
table.x$motif <- paste(
  subseq(table.x$nt.DNAString, 43, 43, NA),
  subseq(table.x$nt.DNAString, 49, 49, NA),
  subseq(table.x$nt.DNAString, 55, 55, NA),
  subseq(table.x$nt.DNAString, 61, 61, NA),
  subseq(table.x$nt.DNAString, 64, 64, NA),
  subseq(table.x$nt.DNAString, 70, 70, NA),
  subseq(table.x$nt.DNAString, 73, 73, NA),
  subseq(table.x$nt.DNAString, 76, 76, NA),
  subseq(table.x$nt.DNAString, 82, 82, NA),
  subseq(table.x$nt.DNAString, 85, 85, NA),
  subseq(table.x$nt.DNAString, 100, 100, NA),
  subseq(table.x$nt.DNAString, 103, 103, NA),
  sep = ""
)

# create a table with total count of motifs 
Total.count<- table.x %>% 
  mutate(Total.count =n()) %>% 
  select(Total.count) %>% 
  unique()

## the sequence below correspond to numbering considering the 4 nt from the index
#table.x <- table.x %>% mutate(Seq = trimws(Seq))
table.x['check']<- NA 
table.x$check <- paste(
  subseq(table.x$nt.DNAString, 14, 14, NA), #T
  subseq(table.x$nt.DNAString, 28, 28, NA), #C
  subseq(table.x$nt.DNAString, 53, 53, NA), #G
  subseq(table.x$nt.DNAString, 58, 58, NA), #C
  subseq(table.x$nt.DNAString, 67, 67, NA), #C
  subseq(table.x$nt.DNAString, 75, 75, NA), #C
  subseq(table.x$nt.DNAString, 94, 94, NA), #C
  subseq(table.x$nt.DNAString, 95, 95, NA), #T
  subseq(table.x$nt.DNAString, 108, 108, NA), #A
  sep = ""
)

## Convert nt.DNAString into character
table.x$nt.DNAString<- as.character(table.x$nt.DNAString)
# 
table.x<- table.x %>% filter(check == "TCGCCCCTA")

## Export clean table with clean barcodes
table.y<- table.x %>% 
  select(Sample, Index, nt.DNAString, motif) ######

# Add number of reads filtered out by parsing conserved sites
Total.count[c("Total.count.2")]<- 
  table.x %>% 
  select(Sample) %>% 
  mutate(n()) %>% 
  select("n()") %>% 
  unique() 

#Total.count <-
#  Total.count %>% 
#  mutate(Discarded = Total.count-Total.count.2)

#Calculate frequency and filter
multi <- table.y %>%
  group_by(motif) %>%
  summarise(n = n()) %>% 
  unique() %>% 
  filter(n >= 3) %>% ### 
  mutate (freq = n / sum(n)) %>% 
  mutate (Index = unique(table.x$Index))

##### Filter table.y based on read count cut off in "multi"

filtered_final<-
  table.y %>% 
  left_join(multi, by = c("motif", "Index")) %>% 
  filter(n >= 3) %>% ## this could also be na.omit
  select("Sample", "Index", "nt.DNAString")

write.table(filtered_final, file = "filtered_final.txt", quote = F, sep="\t",
            row.names = FALSE)

##### Add total count after filtration by read count

Total.count[c("Total.count.3")]<- 
  filtered_final %>% 
  select(Sample) %>% 
  mutate(n()) %>% 
  select("n()") %>% 
  unique() 

write.table(Total.count, file = "Total.count.txt", quote = F, sep="\t",
            row.names = FALSE)

########## Shannon
Shannon <- diversity(multi$n, index="shannon", MARGIN=1, base=exp(1))

write.table(Shannon, 
            file = "Shannon", 
            quote = F, sep="\t",
            col.names = F, 
            row.names = F)

# Export table with barcode frequencies
write.table(multi, 
            file = "multi", 
            quote = F, sep="\t",
            row.names = FALSE)
# For the two write.table above, we removed the txt extension 

#Plot
#makes more colors out of a palette with a limited number
nb.cols <- 4635
mycolors <- colorRampPalette(brewer.pal(8,"Set1"))(nb.cols)

#randomizes the order of colors in the palette
cols = rainbow(4635, s=.5, v=.9)[sample(1:4635,4635)]

multi.plot <- 
  ggplot(multi, aes(x=as.character(Index),y=freq)) + 
  geom_col(aes(fill=motif), show.legend=FALSE, width=0.9) + 
  ylab("Frequency") +
  xlab("Index") +
  scale_fill_manual(values=cols) +
  theme(text=element_text(size=12,face="bold"),
        legend.position = "none",
        plot.title = element_text(size=rel(1.5)),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=rel(1),color = "black"),
        axis.text.y = element_text(size=rel(1),color = "black"),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_line(size=0.5, color = "black"),
        axis.ticks.y = element_line(size = 0.5, color = "black"),
        axis.title.y = element_text(size=rel(1),color = "black"),
        axis.title.x = element_text(size=rel(1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(multi.plot)
ggsave("multi.plot.pdf", width=2,dpi=300)
