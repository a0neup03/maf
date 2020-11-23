---
title: "R Notebook"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

```{r}
require(Biostrings) 
library(BSgenome.Hsapiens.UCSC.hg38)
length(Hsapiens)
#seqnames(Hsapiens)
genome<-BSgenome.Hsapiens.UCSC.hg38
g4_maf<-as.data.frame(filt_maf_df)

g4_maf$snp_left<-abs(g4_maf$start-g4_maf$mcols.loc_start)
g4_maf$snp_right<-g4_maf$end-g4_maf$mcols.loc_end

seq_g4_1 <- BSgenome::getSeq(genome, names = g4_maf$seqnames,
                               start = (g4_maf$start+1),
                               end = g4_maf$start+g4_maf$snp_left-1)

seq_g4_2 <- getSeq(genome, names = g4_maf$seqnames,
                               start = (g4_maf$mcols.loc_start+1),
                               end = g4_maf$mcols.loc_start+g4_maf$snp_right+1)
  
  seq_g4_full <- getSeq(genome, names = g4_maf$seqnames,
                               start = (g4_maf$start+1),
                               end = g4_maf$end+1)
  
  
  head(as.data.frame(seq_g4_1))
head(as.data.frame(seq_g4_2))
head(as.data.frame(seq_g4_full))
g4_seq<-(seq_g4_full)
g4_maf <- g4_maf %>% mutate(G4_id = row_number())
g4_maf_s<-g4_maf[,c(1,2,3,9,14,15,16,17,18,19,21,22,43,44,45,126,127,128,129,130)]
seq_left<-as.data.frame(seq_g4_1)
colnames(seq_left)<-"leftseq"
#last letter of the left seq is the position in question
seq_right<-as.data.frame(seq_g4_2)

rm(seq_g4_1,seq_g4_2)
colnames(seq_right)<-"rightseq"
seq_g4_full<-as.data.frame(seq_g4_full)
colnames(seq_g4_full)<-"seq_g4_full"

seq<-cbind(seq_left,g4_maf_s,as.data.frame(seq_right),seq_g4_full)

seq$G4_ref<-paste0(seq$leftseq,seq$mcols.Reference_Allele,seq$rightseq)
seq$G4_alt<-paste0(seq$leftseq,seq$mcols.Tumor_Seq_Allele2,seq$rightseq)

all_1<-seq

identical(seq$G4_ref,seq$seq_g4_full)



```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.
