library(dplyr)
library(stringr)
library(gplots)
library(ggplot2)
library(plotrix)
library(pheatmap)

library(RColorBrewer)
library(viridis)

variants <- read.table('commonVariants.tsv')
first_line <- readLines("commonVariants.tsv")
first_line <- strsplit(first_line,"\t")[[1]]
colnames(variants) <- first_line

variants$chr_pos_ref_alt <- paste(variants$`#CHROM`,variants$POS,variants$REF,variants$ALT,sep="_")
variants[,c("INFO","ID","QUAL","FILTER","FORMAT","#CHROM","POS","REF","ALT")] <- NULL
rownames(variants) <- variants$chr_pos_ref_alt
variants$chr_pos_ref_alt <- NULL
variants$N_CHR<- NULL

variants[variants==".:."]<-510

for(i in 1:ncol(variants)){
  variants[,i] <-gsub("1:","",variants[,i])
  variants[,i] <- gsub("\\d*,","",variants[,i])
  variants[,i] <- as.numeric(variants[,i])
}

names_df <- data.frame("Condition"=matrix(ncol = ncol(variants), nrow = 1))
colnames(names_df) <- colnames(variants)
for(i in colnames(variants)){
  if(substring(i, 1, 1)=="A"){
    names_df[,i] <- "Atopic dermatitis (A)"
    
  }else if(substring(i, 1, 1) == "N"){
    names_df[,i] <- "Healthy (H)"
  }
  else{
    names_df[,i] <- "Treatment (T)"
  }
}

haplogroup_df <- data.frame("Haplogroup_(HV1)"=matrix(ncol = ncol(variants), nrow = 1))
colnames(haplogroup_df) <- colnames(variants)
for(i in colnames(variants)){
  if(i %in% c("T5","N15","T9")){
    haplogroup_df[,i] <- "A64"
  }
  else if(i %in% c("A7","A4","T6")){
    haplogroup_df[,i] <- "A18"
  }
  else if(i %in% c("A1","N3","N5","N8","T2","N11")){
   haplogroup_df[,i] <- "A68"
  }
  else if(i %in% c("N20","N13","T7")){
    haplogroup_df[,i] <- "A69"
  }
  else if(i %in% c("A8","A6","N6","T1","T4")){
    haplogroup_df[,i] <- "A11"
  }
  else if(i %in% c("N16","A10")){
    haplogroup_df[,i] <- "New haplotype A"
  }
  else if(i %in% c("N14")){
    haplogroup_df[,i] <- "A39"
  }
  else if(i %in% c("N19")){
    haplogroup_df[,i] <- "B1"
  }
  else if(i %in% c("A3","N12","N18")){
    haplogroup_df[,i] <- "New haplotype C"
  }
  else{
    haplogroup_df[,i] <- "C17"
  }
}

mhaplogroup_df <- data.frame("Major_haplogroup"=matrix(ncol = ncol(variants), nrow = 1))
colnames(mhaplogroup_df) <- colnames(variants)
for(i in colnames(variants)){
  if(i %in% c("T5","N15","T9","A7","A4","T6","A1","N3","N5","N8","T2","N11","N20","N13","T7","A8","A6","N6","T1","T4","N16","A10","N14")){
    mhaplogroup_df[,i] <- "A"
    
  }else if(i == "N19"){
    mhaplogroup_df[,i] <- "B"
  }
  else{
    mhaplogroup_df[,i] <- "C"
  }
}

dim(variants)

df <- data.frame("Condition"=t(names_df),"Haplogroup (HV1)"=t(haplogroup_df),"Major_haplogroup"=t(mhaplogroup_df))
dog_names <- colnames(variants)
variants <- as.matrix(variants)

my.breaks <- c(seq(0,255, by=25)) 
my.colors <- c(colorRampPalette(colors = c("black","grey","white"))(length(my.breaks)))

# heatmap with dendogramm for variants and dogs

pdf("heatmap_clusterrow.pdf",         # File name
    width = 8, height = 18, # Width and height in inches
    bg = "white")          # Background color)     
out <- pheatmap(variants,
         cluster_rows =TRUE,
         legend=TRUE,cellwidth = 7,
         cellheight = 4.5,
         fontsize_row = 5,
         fontsize_col = 8,
         fontsize = 7,
         annotation_col =df,
         color=my.colors,
         breaks=my.breaks
         )
dev.off()
