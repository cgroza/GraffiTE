setwd("Google Drive/My Drive/Wheeler_Lab/2024/01-09-GraffiTE_revisions/CORN/Bz_bandage_long_reads/")
table<-read.csv("bz_bandage_family.csv", h = T)
table$Color<-as.factor(table$Color)
table$tclass<-as.factor(table$tclass)
tclass<-levels(table$tclass)
Color<-c("#800000", # DNA/DTA
          "#cc0000", # DNA/DTC
          "#ff1a1a", # DNA/DTH
          "#ff6666", # DNA/DTM
          "#ffb3b3", # DNA/DTT
          "#FFCA00", # DNA/Helitron
          "#87e184", # LTR/Copia
          "#145214", # LTR/Gypsy
          "#00ffbf", # LTR/Unknown
          "#803300", # MITE/DTA
          "#cc5200", # MITE/DTC
          "#ff751a", # MITE/DTH
          "#ffa366", # MITE/DTM
          "#ffd1b3") # MITE/DTT 

index<-as.data.frame(cbind(tclass, Color))
table2<-table[,c("Name", "tclass", "tname")]
table_finish<-merge(table2, index, by = "tclass")
table_finish<-table_finish[,c(2,1,3,4)]
table_finish<-table_finish[order(table_finish$Name),]
write.csv(x = table_finish, file = "bz_custom_cols.csv", quote = F, row.names = F)
