library(dplyr)
# Import file
file1 <-list.files(path = ".", pattern = "\\_multi_final.txt")
# convert to table
table.x<-read.table(file1, header = T)

# Count total number of barcodes 
Richness<-table.x %>% 
  count("motif") %>% 
  select("n")
  
write.table(Richness, 
            file = "Richness.txt", 
            quote = F, sep="\t",
            col.names = F, row.names = F)