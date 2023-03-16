# R Shiny Application for integrative multi-omics bioinformatics
#
# Tailored for reading Galaxy-pipeline outputs
#
# Aim: Provide interactive exploration of the data
######################################################################


#Date: 21.09.2021
#Author: Magnus Arntzen

library(seqinr)

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


getCAZyDetails <- function(cazyTerm) {
  
  #Class
  class <- 'NA'
  if (grepl('CBM', cazyTerm)) class <-'CBM'
  if (grepl('GT', cazyTerm)) class <-'GT'
  if (grepl('CE', cazyTerm)) class <-'CE'
  if (grepl('AA', cazyTerm)) class <-'AA'
  if (grepl('GH', cazyTerm)) class <-'GH'
  if (grepl('PL', cazyTerm)) class <-'PL'
  if (grepl('cohesin', cazyTerm) || grepl('dockerin', cazyTerm) || grepl('SLH', cazyTerm)) class <-'Cellulosome'
  
  #Class number
  classnum <- 'NA'
  if (grepl('CBM', cazyTerm)) classnum <- substring(cazyTerm, 4, 10)
  if (grepl('GT', cazyTerm) || grepl('CE', cazyTerm) || grepl('AA', cazyTerm) || grepl('GH', cazyTerm) || grepl('PL', cazyTerm)) {
    classnum <- substring(cazyTerm, 3, 10)
  }
  if (grepl('_', classnum)) classnum <- unlist(strsplit(classnum, '_'))[1] #remove subclasses
  if (grepl('cohesin', cazyTerm)) classnum <- 'cohesin'
  if (grepl('dockerin', cazyTerm)) classnum <-'dockerin'
  if (grepl('SLH', cazyTerm)) classnum <- 'SLH'
  
  #Pretty text
  prettytext <- 'NA'
  if (class == 'CBM') prettytext <- "CBM - Carbohydrate binding domain"
  if (class == 'GT') prettytext <- "GT - Glycosyl transferases"
  if (class == 'CE') prettytext <- "CE - Carbohydrate esterases"
  if (class == 'AA') prettytext <- "AA - Auxiliary activities"
  if (class == 'GH') prettytext <- "GH - Glycoside hydrolases"
  if (class == 'PL') prettytext <- "PL - Polysaccharide lyases"
  if (class == 'Cellulosome') prettytext <- "Cellulosome"
  
  return(as.data.frame(t(c(CAZyClass = class, CAZyClassNumber = classnum, CAZyText = prettytext))))
}

getKEGG_EC_nr <- function(koTerm, mapping_tbl) {
  
  #koTerm may be semi-colon separated
  
  dat <- data.frame(KO = koTerm) %>% separate_rows('KO', sep=";")
  dat <- dat %>% left_join(mapping_tbl, by= c('KO' = 'KO'))
  
  KO <- paste0(unique(dat$KO), collapse = ";")
  
  #Due to the mapping more redundancies are introduced as one KO can be mapped to several EC and pathways
  EC <- dat %>% 
    select(EC) %>% 
    separate_rows('EC', sep = ';') %>% 
    distinct() %>% 
    rownames_to_column() %>% 
    pivot_wider(names_from = rowname, values_from = EC) %>% 
    unite(col = 'new', sep=";", na.rm = TRUE) %>% 
    deframe()
  #EC <- paste0(EC_distinct, collapse = ";")
  
  KEGG <- dat %>% 
    select(Pathway) %>% 
    separate_rows('Pathway', sep = ';') %>% 
    distinct() %>% 
    rownames_to_column() %>% 
    pivot_wider(names_from = rowname, values_from = Pathway) %>% 
    unite(col = 'new', sep=";", na.rm = TRUE) %>% 
    deframe()

  #KEGG <- paste0(KEGG_distinct, collapse = ";")
    
  return(as.data.frame(t(c(KO = KO, EC = EC, KEGG = KEGG))))
}


get_legend<-function(myggplot){
  #function to extract a legend from a ggplot
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


#JavaScript extensions
jscode <- "
shinyjs.disableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
}

shinyjs.enableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
  tab.unbind('click.tab');
  tab.removeClass('disabled');
}
"

css <- "
.nav li a.disabled {
  background-color: #8da0ab !important;
  color: #333 !important;
  cursor: not-allowed !important;
  border-color: #8da0ab !important;
  font-style: italic;
}"
