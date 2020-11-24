###
##################
library(GenomicRanges)
#library()
library(TCGAbiolinks)
library(data.table)
library(maftools)
library(dplyr)
library(GenomicAlignments)
#select G4 sequences with mutations

#g4_bed<-fread('C:/Users/Aryan Neupane/OneDrive - University of Louisville/Research/Annotate_G4_Clusters_CPP/Annotate_G4_Clusters_CPP/annotation_gencode/',sep= '\t',col.names=c("chromosome","start","end","id","score","strand")) %>%  arrange(chromosome,start,end)%>% distinct()
#make Granges data from the maf file
maf <- GDCquery_Maf("BRCA", pipelines = "muse") %>% read.maf()

g4_bed<-fread('C:/Users/Aryan Neupane/OneDrive - University of Louisville/Research/annotation_G_quad_28feb2020/4_bed.bed',sep= '\t',col.names=c("chromosome","start","end","strand","id")) %>%  arrange(chromosome,start,end)%>% distinct()

##convert maf to granges object

maf_to_granges <- function(maf_df) {
  GenomicRanges::GRanges(
    seqnames = maf_df$Chromosome,
    ranges = IRanges::IRanges(
      start = maf_df$Start_Position,
      end = maf_df$End_Position
    ),
    mcols = maf_df
  )
}
 
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
 

snv_ranges_filter <- function(maf_df,
                              keep_ranges = keep_ranges,
                              bp_window = bp_window)
  {
  # Given a MAF formatted data.frame and a BED regions data.frame; filter out
  # any variants of the MAF df that are not within the BED regions.
  #
  # Args:
  #   maf_df: maf data that has been turned into a data.frame. Can be a maf object
  #           that is subsetted using `@data`.
  #   keep_ranges: BED ranges data.frame with columns: chromosome, start, end
  #             positions in that order or a Genomic Ranges object. If data.frame,
  #             is given, it will be converted to GRanges object
  #   bp_window: how many base pairs away can it be from the BED region to still
  #              be included? Default is 0 bp. This argument gets forwarded
  #              to GenomicRanges::findOverlaps's `maxgap` argument.
  #
  # Returns:
  # The same MAF formatted data.frame with the mutations that lie outside
  # the supplied BED regions filtered out.
  # Turn the MAF sample mutations into a GRanges object
  
  
  #input bed file
  #or
  #g4<-makeGRangesFromDataFrame(g4_bed,keep.extra.columns = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3",ignore.strand = TRUE)
  
  #####get the ratio of G4 in maf file of specific tumors
  
  maf_granges <- maf_to_granges(a)
  
  
  
  if (is.data.frame(keep_ranges)) {
    # Turn the bed regions df into a GRanges object
    keep_ranges <- GenomicRanges::GRanges(
      seqnames = keep_ranges$chromosome,
      ranges = IRanges::IRanges(
        start = keep_ranges$start,
        end = keep_ranges$end
      )
    )
  }
  


  ranges<-subsetByOverlaps(keep_ranges, maf_granges)
  overlap <- findOverlaps(keep_ranges, maf_granges,maxgap = bp_window)
  gr1.matched <- keep_ranges[queryHits(overlap)]
  mcols(gr1.matched) <- cbind.data.frame(
    mcols(gr1.matched),
    mcols(maf_granges[subjectHits(overlap)]));
  
  
  
  g4_maf<-as.data.frame(gr1.matched)
  
  
  # Calculate of ratio of variants in this BED using the @from slot which
  # indicates the indices of the ranges in `maf_ranges` that have overlaps
  # with `keep_ranges`
  
  ratio <- length(overlap@from) / nrow(maf_df)
  
  cat(
    "Ratio of variants in this BED:", ratio, "\n",
    "Ratio of variants being filtered out:", 1 - ratio, "\n"
  )
  
  filt_maf_df <- maf_df[unique(overlap@from), ]
  length_nrow<-length(filt_maf_df)
  print(paste(length_nrow, "Columns intersected"))
  #print(maf_df[subjectHits(overlap)])
  
  
  # Return this matrix with the WXS mutations filtered but WGS the same
  return(g4_maf)
  
  
  
  
}


#input, change start and end position name and save the SNP location as loc_start and loc_end
a<-maf@maf.silent
a$loc_start=a$Start_Position
a$loc_end=a$End_Position
bed_ranges<-keep_ranges(g4_bed)
#g4_SNP<-snv_ranges_filter(a,g4_bed,1)
# Sum up genome sizes

bed_size <- as.numeric(sum(bed_ranges@ranges@width))


#(Tumor_Sample_Barcode ==g4_SNP$Tumor_Sample_Barcode)


# Filter out mutations that are outside of these coding regions.

filt_maf_df <- snv_ranges_filter(a, keep_ranges = g4_bed,-1)

unique(as.data.frame(filt_maf_df$mcols.Tumor_Sample_Barcode))


length(filt_maf_df$)

# Filter to only the sample's mutations
sample_maf_df <- a %>% dplyr::filter(!is.na(a$Tumor_Sample_Barcode)) %>%
  dplyr::filter(Tumor_Sample_Barcode == filt_maf_df$mcols.Tumor_Sample_Barcode)

unique(sample_maf_df$Tumor_Sample_Barcode)



tmb <- filt_maf_df %>%
  dplyr::group_by(
    #TODO: Make this column passing stuff more flexible with some tidyeval maybe
    Tumor_Sample_Barcode = Tumor_Sample_Barcode
  ) %>%
  # Count number of mutations for that sample
  dplyr::summarize(
    mutation_count = dplyr::n(),
    region_size = bed_size,
    tmb = mutation_count / (region_size / 1000000)
  )

unique(sample_maf_df$Tumor_Sample_Barcode)



##filt_maf_df is the file with G4 mutations





