# Script to make master-files containing info from all these Galaxy-pipeline outputs:
#
#  1. List of fasta files, one per bin, from MaxBin2
#  2. Contig coverages, from CoverM
#  3. MAG completeness and contamination, from CheckM
#  4. Gene annotation data, from InterProScan, KoFamScan and dbCAN, one file combined
#  5. Protein quantification, from MaxQuant
#  6. mRNA quantification, from Kallisto
#  7. Taxonomical data, from BAT/CAT
#
#
#  Output:
#  1. List of all contigs, their bins, GC% and coverages
#  2. A tabular file of all genes with bin and contig origin, functional annotation and metaP/T expression values
#
# Date: 08.03.23
# Author: Magnus Arntzen
#
###################################################################################




wd <- getwd()

######    CHANGE THE FOLLOWING TO FIT THE DATA   ######
#------------------------------------------------------

#Required files: (must be present)
contig_abundance_file <- paste0(wd, "/data_from_galaxy/CoverM-Contig.txt")
path_to_bins_fasta <- paste0(wd,"/data_from_galaxy/MaxBin2_ Bins")
annotation_file_name <- paste0(wd, "/data_from_galaxy/InterProScan_KoFamScan_dbCAN.txt")

#Optional files: (Set to NA if not used)
protquant_file_name <- paste0(wd, "/data_from_galaxy/MaxQuant-Protein_groups.txt")
transquant_file_name <- paste0(wd, "/data_from_galaxy/Kallisto-All_abundances_combined.txt")
BAT_file_name <- paste0(wd, "/data_from_galaxy/NEED-TO-RERUN__BAT.txt")
checkM_file_name <- paste0(wd, "/data_from_galaxy/CheckM_lineage_wf.txt")

#Regular expression to extract contig_id from protein_id
rx_contigID_from_proteinID <- "[^_]*_[^_]*"
#Some examples:
#str_extract("NODE_1000_length_967_cov_24.6169_1_884_+", "[^.]*\\.[^_]*")
#str_extract("k141_10012_1359_1805_+", "[^_]*_[^_]*")



###################################################################################

options(stringsAsFactors=F)

library(tidyverse)
library(seqinr)


#################### Read in essential map-data ###################################

KEGG_mapping_file <- paste0(wd, "/../ko00001.keg.mapping.txt")  
ko_mapping_tbl <- read.table(KEGG_mapping_file, header=T, sep="\t", quote="", na.strings = c(""))

ko_mapping_tbl <- ko_mapping_tbl %>% drop_na(Pathway) #Must remove ECs that are only linked to BRITE as these will confuse/introduce errors when trying to map to KEGG-pathways in ViMO

#################### Define some functions ########################################

calculate_contig_GC_percentage <- function(path_to_bins) {
  df<-data.frame(Bin=character(), Contig=character(), GC=double(), Length=integer())
  
  for (f in list.files(path_to_bins, full.names = TRUE)) {
    contig_list <- read.fasta(file = f, seqtype = "DNA")
    bin <- basename(f)
    
    for (n in 1:length(contig_list)) {
      accn <- names(contig_list[n])
      len <- getLength(contig_list[n])
      gc <- GC(contig_list[[n]])
      
      df <- df %>% add_row(Bin = bin, Contig = accn, GC = gc, Length=len)
    }
  }
  return(df)
}


getKEGG_EC_nr <- function(koTerm) {
  #Function to map a K-gene-number to its respective EC and KEGG-pathways
  #Using a mapping table based on the generic downloadable KEGG files
  #Function should return a non-redundant list of ECs and KEGG-Pathways per K-number, and without any NA
  
  
  #koTerm may be list of semi-colon separated K-numbers
  
  dat <- data.frame(KO = koTerm) %>% separate_rows('KO', sep=";")
  dat <- dat %>% left_join(ko_mapping_tbl, by= c('KO' = 'KO'))
  
  KO <- paste0(unique(dat$KO), collapse = ";")
  
  #Due to the mapping, more redundancies are introduced as one KO can be mapped to several EC and pathways
  #We want a non-redundant list, and without NA
  EC <- dat %>% 
    select(EC) %>% 
    separate_rows('EC', sep = ';') %>% 
    distinct() %>% 
    rownames_to_column() %>% 
    pivot_wider(names_from = rowname, values_from = EC) %>% 
    unite(col = 'new', sep=";", na.rm = TRUE) %>% 
    deframe()

  KEGG <- dat %>% 
    select(Pathway) %>% 
    separate_rows('Pathway', sep = ';') %>% 
    distinct() %>% 
    rownames_to_column() %>% 
    pivot_wider(names_from = rowname, values_from = Pathway) %>% 
    unite(col = 'new', sep=";", na.rm = TRUE) %>% 
    deframe()
  
  return(as.data.frame(t(c(KO = KO, EC = EC, KEGG = KEGG))))
}






###################################################################################
#####################    Contig-file    ###########################################
###################################################################################




############### Load MAGs and calculate GC content ################################

print("Reading MAGs and calculating GC%")

df.GC <- calculate_contig_GC_percentage(path_to_bins_fasta)
df.GC <- df.GC %>%
  mutate(Bin = gsub('.fasta|.faa|.fna', '', Bin)) %>%
  mutate(GC = GC*100)


##################### Adding contig coverages #####################################

print("Reading contig coverages")


contigs <- read.table(contig_abundance_file, header=T, sep="\t", quote="")
contigs_avg <- contigs %>%
  pivot_longer(cols = -Contig) %>%
  group_by(Contig) %>%
  summarise(mean  = mean(value))
names(contigs_avg) <- c('Contig', 'Coverage')

df.GC <- df.GC %>%
  left_join(contigs_avg, by = c("Contig" = "Contig"))

##################### Adding taxonomy #############################################

if (!is.na(BAT_file_name)) {
  print("Reading lineage from BAT")
  
  taxa <- read.table(BAT_file_name, header=T, sep="\t", quote="", comment.char = "", colClasses = "character")
  cols <- colnames(taxa)
  cols[1] = "Bin"
  names(taxa) <- cols
  taxa <- taxa %>%
    select(-classification, -reason, -lineage, -lineage.scores, -full.lineage.names) %>%
    pivot_longer(cols=-1) %>%
    group_by(Bin) %>%
    summarise(Lineage = paste0(gsub(" \\(.*","",value), collapse = ";"))
    
  df.GC <- df.GC %>%
    left_join(taxa, by = c("Bin" = "Bin"))
  
}

##################### Adding MAG completeness and contamination ####################

print("Reading CheckM data")

if (!is.na(checkM_file_name)) {
  
  checkM <- read.table(checkM_file_name, header=T, sep="\t", quote="", comment.char = "", colClasses = "character")
  cols <- colnames(checkM)
  cols[1] = "Bin"
  names(checkM) <- cols
  checkM <- checkM %>%
    select(Bin, Completeness, Contamination,	Strain.heterogeneity)

  df.GC <- df.GC %>%
    left_join(checkM, by = c("Bin" = "Bin"))
  
}

##################### Writing contig-file ##########################################

bin_contig_coverage_file = "Bin_taxonomy_coverages.txt"
print(paste("Writing overview of bins, contigs and coverages to", bin_contig_coverage_file))

write_delim(df.GC, bin_contig_coverage_file, delim = "\t")










###################################################################################
#####################    Master-file    ###########################################
###################################################################################



############# Load InterProScan annotation file and add Contigs and Bins ##########

print("Reading InterProScan results and annotate with bin information")

tbl <- read.table(annotation_file_name, header=T, sep="\t", quote="")

tbl <- tbl %>%
  mutate(Contig = str_extract(Accn, rx_contigID_from_proteinID)) %>%
  left_join(df.GC %>% select(Contig, Bin), by = c("Contig" = "Contig")) %>%
  mutate(Bin = replace(Bin, is.na(Bin), 'Unbinned'))




############# Add KEGG and EC from KoFamScan-mapping ##############################

#Since 2022, InterProScan does no longer provide KO-pathways and EC annotations,
#hence, the KEGG and EC columns are often empty, which will affect the KEGG analysis
#in this tool. We need to map these back using the ko00001.keg files, given that the data
#contains KO info from KoFamScan


if ('KoFamScan' %in% colnames(tbl)) {
  
  #The App uses 'KO' as name for KoFamScan
  tbl <- tbl %>% rename(KO = KoFamScan)
  
  
  if (nrow(tbl %>% drop_na(KEGG)) == 0) {

    print("Remapping KEGG/EC from KoFamScan. This may take a few minutes")
    
    #dynamic programming: foreach K-number, add ECs and KOs, but non-redundant and without NA.
    tbl <- tbl %>% select(-EC, -KEGG)
    tbl$kannot <- lapply(tbl$KO, function(x) getKEGG_EC_nr(x))
    tbl$EC <- unlist(lapply(tbl$kannot, "[[", 2))
    tbl$KEGG <- unlist(lapply(tbl$kannot, "[[", 3))
    tbl <- tbl %>% select(-kannot)
  }
}





######## Load protein quantification file + filter and annotate
if (!is.na(protquant_file_name)) {
  
  print("Reading MaxQuant proteomics results and annotate the InterProScan results accordingly")
  
  tbl_protquant <- read.table(protquant_file_name, header=T, sep="\t", quote="", na.strings = "NaN")
  
  #Filter MaxQuant output, common base filtering: 1) remove contaminants, reverse hits and site ids, and at least one experiments is non-NaN
  df.protquant <- tbl_protquant %>%
    select(Majority.protein.IDs, Only.identified.by.site, Reverse, Potential.contaminant, contains('LFQ.intensity')) %>%
    rename_with(~ gsub('LFQ.intensity.', 'metaP_LFQ_', .x, fixed = TRUE)) %>%
    filter(!Reverse=='+') %>%
    filter(!Only.identified.by.site=='+') %>%
    filter(!Potential.contaminant=='+') %>%
    select(-Only.identified.by.site, -Reverse, -Potential.contaminant) %>%
    mutate(across(contains('metaP_LFQ_'), log2)) %>%
    na_if("-Inf") %>%
    mutate(across(contains('metaP_LFQ_'), function(x) 2^x)) %>%
    filter(if_any(contains('metaP_LFQ_'), ~ !is.na(.)))
  
  #Open up the semicolon-separated IDs and give the quantitative data to all of them
  df.protquant_expanded <- df.protquant %>%
    separate_rows(Majority.protein.IDs, sep = ";")
  
  #Append this to the InterProScan table
  tbl <- tbl %>%
    left_join(df.protquant_expanded, by = c("Accn" = "Majority.protein.IDs"))
}






######## Load transcript quantification file + annotate
if (!is.na(transquant_file_name)) {
  
  print("Reading Kallisto transscriptomics results and annotate the InterProScan results accordingly")
  
  tbl_transquant <- read.table(transquant_file_name, header=T, sep="\t", quote="", na.strings = "NaN")
  
  #Filter Kallisto files:  at least one experiments is non-NaN
  df.transquant <- tbl_transquant %>%
    select(target_id, contains('_tpm')) %>%
    rename_with(~ paste0("metaT_TPM_", .x), contains('_tpm')) %>%
    rename_with(~ gsub('_tpm', '', .x, fixed = TRUE)) %>%
    rename_with(~ gsub('.fastq', '', .x, fixed = TRUE)) %>%
    mutate(across(contains('_tpm'), log2)) %>%
    na_if("-Inf") %>%
    mutate(across(contains('_tpm'), function(x) 2^x)) %>%
    filter(if_any(contains('_tpm'), ~ !is.na(.)))
  
  #Open up the semicolon-separated IDs (if any) and give the quantitative data to all of them
  df.transquant_expanded <- df.transquant %>%
    separate_rows(target_id, sep = ";")
  
  #Append this to the InterProScan table
  tbl <- tbl %>%
    left_join(df.transquant_expanded, by = c("Accn" = "target_id"))
}




######## Save final file
final_file = "All_genes_functional_quantitative_data.txt"
print(paste("Writing final annotation file of all genes with functional and quantitative data to", final_file))

write_delim(tbl, final_file, delim = "\t")

