


#------------------------------------------
#files required:
#G4 annotation file : bed_annot1.1
#reference region: eg promoters, exons, introns, splicing sites, experimental regions
# functions : 1. intersection and get genes for further enrichment
#enrichment analysis function

# get genes related to different GO terms
# get G4 related to the genes

#here we go

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(dplyr)
library(GenomicAlignments)
library(enrichplot)
library(GOstats)
library(org.Hs.eg.db)
library(enrichplot)
library(GOstats)
library(tidyverse)

library(clusterProfiler)
###

keep_ranges<-function(keep_ranges){
  #if (is.na(flank)) {
  if (is.data.frame(keep_ranges)) {
    # Turn the bed regions df into a GRanges object
    keep_ranges <- GenomicRanges::GRanges(
      seqnames = keep_ranges$chromosome,
      ranges = IRanges::IRanges(
        start = keep_ranges$start,
        end = keep_ranges$end
      ),
      mcols = keep_ranges
    )
    
  }
  
}
#maf_df is the region to compare with, example G4 region present in intron, exon etc
maf_df<- #  make these files
  G4_bed<-bed_annot1.1
snv_ranges_filter <- function(maf_df,G4_bed,
                              bp_window = bp_window)
  
  
  G4_bed_granges<-keep_ranges(G4_bed)
maf_df_granges<-keep_ranges(maf_df)

ranges<-subsetByOverlaps(G4_bed_granges, maf_df_granges)
if (is.na(bp_window)==TRUE) {
  bp_window= 1}
overlap <- findOverlaps(G4_bed_granges, maf_df_granges,minoverlap=bp_window)
gr1.matched <- G4_bed_granges[queryHits(overlap)]
mcols(gr1.matched) <- cbind.data.frame(
  mcols(gr1.matched),
  mcols(maf_df_granges[subjectHits(overlap)]),gr1.matched);

g4_maf<-as.data.frame(gr1.matched)
# Calculate of ratio of variants in this BED using the @from slot which
# indicates the indices of the ranges in `maf_ranges` that have overlaps
# with `keep_ranges`

ratio <- length(overlap@from) / nrow(maf_granges)

cat(
  "Ratio of variants in this BED:", ratio, "\n",
  "Ratio of variants being filtered out:", 1 - ratio, "\n"
)

filt_maf_df <- a[unique(overlap@from), ]
length_nrow<-length(filt_maf_df)
print(paste(length_nrow, "Columns intersected"))
#print(maf_df[subjectHits(overlap)])

#table(overlap@from)
# Return this matrix with the WXS mutations filtered but WGS the same
return(g4_maf)
}

##

# get regions
txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene
genes = transcriptsByOverlaps(ranges = bed, x = db, maxgap=window, columns=c('tx_id', 'gene_id'))

orgdb = org.Hs.eg.db
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
head(bed_annot1.1)




#make regex for any term, or change use multiple words as c("a","b")
regex <- paste0(sprintf("(?=.*%s)", "same.region"), collapse = '')

annot_to_select<-"PROMOTER"
#"PROMOTER"   "Intron"     "Intergenic" "3' UTR"     "Downstream" "EXON"       "5' UTR" 
bed_annot1_selected<-bed_annot1.1[(grepl(regex,bed_annot1.1$rowid,perl=TRUE)),] 

genes_selected<-bed_annot1_selected %>% filter(annot==annot_to_select) %>% dplyr::select(geneId)

universe=unique(bed_annot1_selected$geneId)

#universeGeneIds <- keys(org.Hs.eg.db)
go_enrich <- enrichGO(gene = genes_selected$geneId,
                      universe = universe,
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)


list1<-herum %>% 
  mutate(genes = strsplit(as.character(test), ",")) %>%
  unnest(test) %>% select(test)%>%slice_sample() 
herum<- temp1%>% group_by(V14) %>%select(V13)%>% group_split()
result<-as.data.frame(go_enrich@result)




dft %>%
  rowwise() %>%
  mutate( g = sum(t1 %in% k)) %>%
  filter( g > 0) %>%
  select(t1, cnt)
