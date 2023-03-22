###################################################################################
#
# ViMO (Visualizer for multi-omics)
#
# An R Shiny Application for integrative multi-omics bioinformatics
#
# Inputs: 
#
#  1. Masterfile of all genes with quantitative and functional annotation (tsv)
#  2. File with all contig, their bins, GC% and coverages (tsv)
#
# Aim: Provide interactive exploration of the data
#
#
# Date: 26.09.22
# Author: Magnus Arntzen
#
###################################################################################


#You want to increase this if uploading Masterfile and Contigfile > 200 Mb.
options(shiny.maxRequestSize=200*1024^2)


######################## -  LIBRARIES AND GLOBALS - #####################################
options(stringsAsFactors=F)

library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggthemr)
library(RColorBrewer)
library(dichromat)
library(grid)
library(gridExtra)
library(DT)
library(KEGGREST)
library(shiny)
library(shinythemes)
library(shinyjs)
library(MetQy)
library(pheatmap)
library(magick)
wd <- getwd()

options(dplyr.summarise.inform = FALSE)
options(repos = BiocManager::repositories())


######################## -  EXTERNAL FILES AND GLOBALS- #####################################
source("functions.R")


#No need to change anything here.
#ViMO allows for uploading your own files at run-time

KEGG_generic_file <- paste0(wd, "/ko00001.keg.txt")  
KEGG_mapping_file <- paste0(wd, "/ko00001.keg.mapping.txt")

masterfile <- ''
contigfile <- ''
example_masterfile <- paste0(wd, "/www/All_genes_functional_quantitative_data.txt")
example_contigfile <- paste0(wd, "/www/Bin_taxonomy_coverages.txt")




######################## -  LOADING AND PREPROCESSING - #####################################


# Load KEGG generic data for App functionality
ko_tbl <- read.table(KEGG_generic_file, header=T, sep="\t", quote="")
ko_mapping_tbl <- read.table(KEGG_mapping_file, header=T, sep="\t", quote="", na.strings = c(""))
k_numbers_connected_to_pathways <- ko_mapping_tbl %>%
                            drop_na(Pathway) %>%
                            select(KO) %>%
                            distinct() %>%
                            deframe()


######################## -  APP UI - #####################################

#Define UIs

ui_dataset_details <- fluidRow(
  wellPanel(
    h3("Dataset overview"),
    fileInput("masterfile", "Choose Masterfile", accept = ".txt"),
    fileInput("contigfile", "Choose Contigfile", accept = ".txt"),
    p("See the help section for an explanation of the file formats."),
    p("Current file-size upload limits are set to 200 Mb."),
    actionButton("loadButton", "Load files", class = "btn-success"),
    actionButton("loadButton_example", "Use example files", class = "btn-success"),
    
    hr(),
    
    htmlOutput("filesSummary"),
    
    #actionButton("STOPButton", "Close App", icon = icon("fa-solid fa-octagon-exclamation")),
    
  )
)



ui_cazy_functional_analysis <- fluidPage(width = 12,
  wellPanel(
    tags$style(type='text/css', ".selectize-input { padding: 2px; min-height: 0;} .selectize-dropdown { line-height: 10px; }"),
    selectInput(inputId = "cazyGroup", label = "Choose a CAZy functional group", choices = c("NA"))
  ),
  fluidRow(
    h3(textOutput("cazyCaption", container = span)),
    fluidRow(
      column(4, offset = 4,
        selectInput(inputId = "sample_MT", label = "MetaT sample", choices = c("NA")),
      ),
      column(4,
        selectInput(inputId = "sample_MP", label = "MetaP sample", choices = c("NA")),
      )
    ),
    fluidRow(
      column(4, offset = 8,
        selectInput(inputId = "yaxis_MP", label = "y-axis", choices = c("Summed LFQ", "Protein count")),
      )
    ),
    fluidRow(
      column(4,
        h4("Meta-genomics"),
        plotOutput(outputId = "cazyPlot_MG"),
        downloadButton("download_cazyPlot_MG", label = "PDF")
      ),
      column(4,
        h4("Meta-transcriptomics"),
        plotOutput(outputId = "cazyPlot_MT"),
        downloadButton("download_cazyPlot_MT", label = "PDF")
      ),
      column(4,
        h4("Meta-proteomics"),
        plotOutput(outputId = "cazyPlot_MP"),
        downloadButton("download_cazyPlot_MP", label = "PDF")
      )
    ),
    fluidRow(
      column(12, 
        plotOutput(outputId = "cazyPlot_legend", inline = TRUE),# height = "100px"),
        hr(),
        h5("The barplots show the number/abundance of genes/transcripts/proteins per bin and per sample."),
        h5("The heatmaps below show the abundance of transcripts/proteins per bin but for all samples simultaneously to highlight e.g. changes over time.")
      )
    ),
    hr(),
    fluidRow(
      column(12, 
        h3("Functional abundance per sample per MAG"),
        selectInput(inputId = "cazyPlot_heatmap_binSelector", label = "MAG:", choices = c("NA")),
        plotOutput(outputId = "cazyPlot_heatmap_MT", inline = TRUE),
        downloadButton("download_cazyPlot_MT_heatmap", label = "PDF"),
        hr(),
        plotOutput(outputId = "cazyPlot_heatmap_MP", inline = TRUE),
        downloadButton("download_cazyPlot_MP_heatmap", label = "PDF")
      )
    )
  )
)

ui_kegg_functional_analysis <- fluidPage(width = 12,
   wellPanel(
     selectInput(inputId = "keggGroup", label = "Choose a KEGG pathway", choices = c("NA"))
   ),
  fluidRow(
    h3(textOutput("keggCaption", container = span)),
    fluidRow(
      column(4, offset = 4,
             selectInput(inputId = "kegg_sample_MT", label = "MetaT sample", choices = c("NA")),
      ),
      column(4,
             selectInput(inputId = "kegg_sample_MP", label = "MetaP sample", choices = c("NA")),
      )
    ),
    fluidRow(
      column(4, offset = 8,
             selectInput(inputId = "kegg_yaxis_MP", label = "y-axis", choices = c("Summed LFQ", "Protein count")),
      )
    )
  ),
  fluidRow(
    column(4,
       h4("Meta-genomics"),
       plotOutput(outputId = "keggPlot_MG"),
       downloadButton("download_keggPlot_MG", label = "PDF")
    ),
    column(4,
       h4("Meta-transcriptomics"),
       plotOutput(outputId = "keggPlot_MT"),
       downloadButton("download_keggPlot_MT", label = "PDF")
    ),
    column(4,
       h4("Meta-proteomics"),
       plotOutput(outputId = "keggPlot_MP"),
       downloadButton("download_keggPlot_MP", label = "PDF")
    )
  ),
  fluidRow(
    plotOutput(outputId = "keggPlot_legend", inline = TRUE)#, height = "100px")
  ),
  fluidRow(
      div(DT::dataTableOutput("keggTable"), style ="font-size: 90%; height: 75%"),
      em("Click a pathway to download and highlight the corresponding KEGG pathway map"),
      hr()
  ),
  fluidRow(
    column(7,
      h4(htmlOutput("ERRORtext_KEGGmap")),
      imageOutput(outputId = "keggPathway", inline = TRUE), #height = "600px"),
      hr(),
      em("NB: Note that the scale and colors used in the KEGG map can deviate from the heatmaps to the right as this depend on color swath chosen in settings and the number of bins in the data")
    ),
    column(5,
      checkboxInput(inputId = "keggplot_check_download_kegg", label = "Download and display KEGG pathway", value = TRUE),
      radioButtons(inputId = "keggplot_annotation_radio", label = "Mark KEGG pathway with", choiceNames =  c("Meta-genomics", "Meta-transcriptomics (average abundance across samples)", "Meta-proteomics (average abundance across samples)"), choiceValues =  c("MetaG", "MetaT", "MetaP")),
      selectInput(inputId = "keggPlot_heatmap_binSelector", label = "MAG:", choices = c("NA")),
      plotOutput(outputId = "keggPlot_heatmap_MT", inline = TRUE),
      downloadButton("download_keggPlot_MT_heatmap", label = "PDF"),
      hr(),
      plotOutput(outputId = "keggPlot_heatmap_MP", inline = TRUE),
      downloadButton("download_keggPlot_MP_heatmap", label = "PDF")
    )
  )
)

ui_MAGs <- fluidPage(width = 12,
  fluidRow(
    h3("Metagenome-assembled genomes (MAGs)"),
    h4(textOutput("MAGplot_info", container = span)),
    fluidRow(
      column(2,
        sliderInput(inputId = "MAG_slider_completeness", label = "Filter to MAGs with Completeness above", min = 0, max = 100, value = 90, step = 1, post="%", ticks = FALSE),
      ),
      column(2,
        sliderInput(inputId = "MAG_slider_contamination", label = "Filter to MAGs with Contamination below", min = 0, max = 100, value = 5, step = 1, post="%", ticks = FALSE),
      ),
      column(8)
    ),
    plotOutput(outputId = "MAGPlot", height = "600px"),
    downloadButton("downloadMAGPlot", label = "PDF"),
    h5("The plot shows the GC% and coverage of all contigs in the metagenomics data. The size reflects the contig length and the color shows which MAG/bin the contig was assigned to."),
    h5("The plot will never show more then 50 MAGs for visualization purposes. If the table below has >50 MAGs the plot will be limited to the 50 most abundant ones."),
    br(),
    br(),
    div(DT::dataTableOutput("MAGTable"), style ="font-size: 90%; height: 75%"),
    em("Click a bin to highlight it in the plot"),
    hr(),
    h3("MAG abundances"),
    plotOutput(outputId = "MAGabundancePlot_heatmap_MT", inline = TRUE),
    downloadButton("download_MAGabundancePlot_MT_heatmap", label = "PDF"),
    hr(),
    plotOutput(outputId = "MAGabundancePlot_heatmap_MP", inline = TRUE),
    downloadButton("download_MAGabundancePlot_MP_heatmap", label = "PDF")
  )
)

ui_MCF <- fluidPage(width = 12,
    fluidRow(
      column(1,
        actionButton("calculate_mcf_button", "Calculate mcf", class = "btn-success")
      ),
      column(11,
        em("NB: Calculation of Module completion fraction (mcf) is time-consuming and more so for large datasets with many bins. For the example files (9 bins, 15.656 genes), this will take 3-4 minutes and ViMO will be unresponsive while calculating. The plot will be updated when the process is finished.")
      )
    ),
    br(),
    br(),
    wellPanel(
      selectInput(inputId = "mcf_keggGroup", label = "Filter by a KEGG pathway", choices = c("NA")),
      numericInput(inputId = "mcf_module_cutoff", label = "Only show modules where at least one mcf value is at/above", value = 0.5, min = 0.0, max = 1.0, step = 0.1)
    ),
    fluidRow(
     plotOutput(outputId = "mcfPlot", height = "1200px"),
     h5("The plot shows the Module completion fraction (mcf) for all bins. This is indicative of the metabolic potential of each bin/MAG"),
     column(1,
      downloadButton("downloadMCFPlot", label = "PDF")
     ),
     column(1,
      downloadButton("downloadMCFPlot_data", label = "tsv")
     )
   )
)

ui_Settings <- fluidPage(width = 12,
   wellPanel(
     h3("Settings"),
     selectInput(inputId = "gg_theme", label = "Theme for ggplots", choices = c('flat', 'flat dark', 'dust', 'light', 'earth', 'fresh', 'chalk', 'lilac', 'carrot', 'pale', 'copper', 'grape', 'greyscale', 'sky', 'solarized', 'grass', 'sea', 'camouflage'), selected = 'flat'),
     textOutput("gg_themeCaption", container = span)
     #shinythemes::themeSelector()
   )
)

ui_Help <- fluidPage(width = 12,
                     fluidRow(
                       includeHTML(paste0(wd, "/helptext.html"))
                     )
)

ui_main <- navbarPage(title=div(img(src="ViMO_logo.png", width="150")), theme = shinytheme("cerulean"), windowTitle = "ViMO",
    tags$head(
      tags$style(
        HTML("#shiny-notification-panel {
              top: 0;
              bottom: unset;
              left: 0;
              right: 0;
              margin-left: auto;
              margin-right: auto;
              width: 100%;
              max-width: 450px;
            }
            .navbar {
              min-height:60px !important;
            }
         ")
      )
    ),
                      
    useShinyjs(),
    extendShinyjs(text = jscode, functions = c("disableTab", "enableTab")),
    inlineCSS(css),
                      
    tabPanel("Load data", ui_dataset_details),
    tabPanel("MAGs", ui_MAGs),
    tabPanel("CAZy", ui_cazy_functional_analysis),
    tabPanel("KEGG", ui_kegg_functional_analysis),
    tabPanel("KEGG-Modules", ui_MCF),
    tabPanel("Settings", ui_Settings),
    tabPanel("Help", ui_Help)
)





######################## -  APP SERVER - #####################################


#Define server
server <- function(input, output, session) {

  #Server config
  session$onSessionEnded(function() {
    stopApp()
  })
  
  
  #Some start-up settings
  js$disableTab("MAGs")
  js$disableTab("CAZy")
  js$disableTab("KEGG")
  js$disableTab("KEGG-Modules")
  
  
  #Colors
  global_swatch = NULL
  
  #Runtime
  mcf_running = FALSE
  
  
  #App functionality, will be adjusted based on loaded files
  do_metaP = TRUE
  do_metaT = TRUE
  do_lineage = TRUE
  do_checkM = TRUE
  
  
  #The main data object is reactive. When changed, it will trigger changes in everything
  data <- reactiveValues(
    masterfile = NULL,
    contigfile = NULL,
    master_tbl = NULL,
    contig_tbl = NULL,
    cazymes = NULL,
    cazymes_expanded = NULL,
    kegg_proteins = NULL,
    kegg_proteins_expanded = NULL
    )
  
  
  
  
  
  #Create all the needed tables for the app once "data" changes (typically when loaded)
  updateTables <- function() {
    print("update tables")
    
    
    #Sanity check of the files
    required_columns_masterfile <- c('Accn', 'Bin', 'KEGG', 'KO', 'CAZy') %in% colnames(data$master_tbl)
    required_columns_contigfile <- c('Contig', 'Bin', 'Length', 'Coverage', 'GC') %in% colnames(data$contig_tbl)
    if (FALSE %in% required_columns_masterfile) {
      showNotification("Masterfile missing required columns, see help section", type = "error", duration = NULL)
      clearAllData()
      return()
    }
    if (FALSE %in% required_columns_contigfile) {
      showNotification("Contigfile missing required columns, see help section", type = "error", duration = NULL)
      clearAllData()
      return()
    }
    do_metaP <<- length(colnames(data$master_tbl %>% select(contains("metaP_LFQ")))) > 0
    do_metaT <<- length(colnames(data$master_tbl %>% select(contains("metaT_TPM")))) > 0
    do_lineage <<- 'Lineage' %in% colnames(data$contig_tbl)
    do_checkM <<- !FALSE %in% (c('Completeness', 'Contamination', 'Strain.heterogeneity') %in% colnames(data$contig_tbl))
    
    print(paste("MetaP:", do_metaP))
    print(paste("MetaT:", do_metaT))
    print(paste("Lineage:", do_lineage))
    print(paste("CheckM:", do_checkM))
    
    
    
    #Since 2022, InterProScan does no longer provide KO-pathways and EC annotations,
    #hence, the KEGG and EC columns are often empty, which will affect the KEGG analysis
    #in this tool.
    if (nrow(data$master_tbl %>% drop_na(KEGG)) == 0) {
      showNotification("No KEGG-Pathway data available from InterProScan. See the help section and use the script makeMasterFile.R to re-map KEGG and EC data from KoFamScan-data", type = "error", duration = NULL)
      clearAllData()
      return()
    }


    
    # Using the Masterfile-table and Contigfile-table, we will extract some sub-tables for easier analysis:
    # The purpose of having these sub-tables is to avoid recalculating them many times - saves time
    #
    # 1. cazymes
    # 2. cazymes_expanded
    # 3. kegg_proteins
    # 4. kegg_proteins_expanded
    
    
    ### Cazymes - a list of all cazymes and which bin they belong to + quantitative data
    data$cazymes <- data$master_tbl %>%
      select(Bin, CAZy, contains("metaT_TPM"), contains("metaP_LFQ")) %>%
      filter(CAZy != "")
      
    ### Cazymes_expanded - open the semicolon separation and annotate more details
    data$cazymes_expanded <- data$cazymes %>% separate_rows(CAZy, sep = ";")
    data$cazymes_expanded$annot <- lapply(data$cazymes_expanded$CAZy, function(x) getCAZyDetails(x))
    data$cazymes_expanded$CAZyClass <- unlist(lapply(data$cazymes_expanded$annot, "[[", 1))
    data$cazymes_expanded$CAZyClassNumber <- unlist(lapply(data$cazymes_expanded$annot, "[[", 2))
    data$cazymes_expanded$CAZyText <- unlist(lapply(data$cazymes_expanded$annot, "[[", 3))
    data$cazymes_expanded <- data$cazymes_expanded %>% select(-annot)

    
    
    ### kegg_proteins - a list of all proteins with KEGG annotation and which bin they belong to + quantitative data
    data$kegg_proteins <- data$master_tbl %>%
      select(Bin, KEGG, EC, KO, contains("metaT_TPM"), contains("metaP_LFQ")) %>%
      filter(KEGG != "")
      
    ### kegg_proteins_expanded - open the semicolon separation and annotate with more detailed KEGG info
    data$kegg_proteins_expanded <- data$kegg_proteins %>% separate_rows(KEGG, sep = ";")
    data$kegg_proteins_expanded <- data$kegg_proteins_expanded %>%
      left_join(ko_tbl %>% select(ID, Name, Parent.term), by = c("KEGG" = "ID"))
    

    #Some extra columns for counting purposes
    data$master_tbl <- data$master_tbl %>%
      mutate(metaP_expressed = if_any(contains('metaP_LFQ'), ~ !is.na(.))) %>%
      mutate(metaT_expressed = if_any(contains('metaT_TPM'), ~ !is.na(.)))

    #Update color swatches
    setThemeAndCalculateColors('flat')
    
    
    #Update UI based on data
    #CAZy-page
    updateSelectInput(session, "cazyGroup", choices = unique(data$cazymes_expanded$CAZyText))
    updateSelectInput(session, "sample_MT", choices = unique(colnames(data$master_tbl %>% select(contains("metaT_TPM")))))
    updateSelectInput(session, "sample_MP", choices = unique(colnames(data$master_tbl %>% select(contains("metaP_LFQ")))))
    updateSelectInput(session, "cazyPlot_heatmap_binSelector", choices = unique(data$contig_tbl$Bin))
    #KEGG-page
    updateSelectInput(session, "keggGroup", choices = unique(data$kegg_proteins_expanded$Parent.term))
    updateSelectInput(session, "kegg_sample_MT", choices = unique(colnames(data$master_tbl %>% select(contains("metaT_TPM")))))
    updateSelectInput(session, "kegg_sample_MP", choices = unique(colnames(data$master_tbl %>% select(contains("metaP_LFQ")))))
    updateSelectInput(session, "keggPlot_heatmap_binSelector", choices = unique(data$contig_tbl$Bin))
    #MCF-page
    #This will be updated after mcf calculation due to different content for modules
    #MAG-Page
    if (do_checkM)
      updateSliderInput(session, "MAG_slider_contamination", max = ceiling(max(data$contig_tbl$Contamination)))
    

    #Enable analysis
    js$enableTab("MAGs")
    js$enableTab("CAZy")
    js$enableTab("KEGG")
    if ('KO' %in% colnames(data$master_tbl))
      js$enableTab("KEGG-Modules")
    else
      showNotification("No 'KO'-column detected. KEGG-Modules analysis disabled", type = "warning", duration = NULL)
    
  }
  

  
  
  
  
    
  
  #After data is loaded and tables are updated, we need to calculate a color swatch suitable for the
  #number of bins in the data set.
  #
  # NB: Note the <<- to change the parent variable 'global_swatch'
  #
  setThemeAndCalculateColors <- function(theme_string) {
    print("recalulate color swatch")

    #Set base theme. This can be tweaked by some UI options if wanted
    ggthemr(theme_string, type="outer", layout="scientific", spacing=2, text_size = 14)
    n_cols <- length(unique(data$contig_tbl$Bin))
    
    #Provides 8 usable bin-colors. Nine in swatch but first one is special
    global_swatch <<- as.vector(swatch())
    
    if (n_cols > 8 & n_cols < 20) {
      #Provides 19 + 1 special
      global_swatch <<- c(global_swatch, brewer.pal(11, "BrBG"))
    }

    if (n_cols > 19 & n_cols < 71) {
      #Provides 70 + 1 special
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      global_swatch <<- c("#555555", unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))))
    }

    if (n_cols > 70 & n_cols < 503) {  #This is becoming a headache!
      #Provides 502 + 1 special
      global_swatch <<- c("#555555", colors(distinct = TRUE))
    }

    if (n_cols > 502) { #AAAhh!...What to do here? There may be datasets with >1500 bins
      #Not possible to generate this many unique colors
      #Not possible to use ggplot with not-unique colors
      #Not possible to draw the plot at all 
      #Not possible to interpret the plot anyhow
      
      #Can we filter to the most abundant MAGs?
      
      #Can we say that ViMO do not support graphing this many MAGs? The users are referred to manual analysis in R.
      
      print(">502 Bins in the data. R, ggplot and ViMO have difficultis with that any. R will likely crash...")
    }


    names(global_swatch) <<- c(NA, sort(unique(data$contig_tbl$Bin)))
    
    set_swatch(global_swatch)

    print(paste0("New theme set:", theme_string, ", Color count:", length(global_swatch)))
  }

  
  
    
  clearAllData <- function() {
    print("Clearing all data...")
    data$masterfile = NULL
    data$contigfile = NULL
    data$master_tbl = NULL
    data$contig_tbl = NULL
    data$cazymes = NULL
    data$cazymes_expanded = NULL
    data$kegg_proteins = NULL
    data$kegg_proteins_expanded = NULL
    
    global_swatch <<- NULL
    do_metaP <<- TRUE
    do_metaT <<- TRUE
    do_lineage <<- TRUE
    do_checkM <<- TRUE
    
    
    js$disableTab("MAGs")
    js$disableTab("CAZy")
    js$disableTab("KEGG")
    js$disableTab("KEGG-Modules")
    
  }
  
  
  # observeEvent(input$STOPButton, {
  #   stopApp()
  # })
  

  #The load buttons trigger a change in the main data object
  observeEvent(input$loadButton, {
    clearAllData()
    
    print("load_button")
    print(paste("masterfile:", input$masterfile$datapath))
    print(paste("contigfile:", input$contigfile$datapath))
    
    #check that the two files are set. Otherwise pop up message to user
    if (is.null(input$masterfile$datapath)) {
      showNotification("No Masterfile set", type = "error")
      return()
    }
    if (is.null(input$contigfile$datapath)) {
      showNotification("No Contigfile set", type = "error")
      print("no contig")
      return()
    }
    
    showNotification("Loading...", type = "message", duration = 2)
    
    data$master_tbl <- read.table(input$masterfile$datapath, header=T, sep="\t", quote="", colClasses=c(Bin="character"))
    data$contig_tbl <- read.table(input$contigfile$datapath, header=T, sep="\t", quote="", colClasses=c(Bin="character"))
    data$masterfile <- paste(input$masterfile$name, '(uploaded)')
    data$contigfile <- paste(input$contigfile$name, '(uploaded)')
    
    
    showNotification("Files loaded. Calculating subsets...", type = "message", duration = 5)
    updateTables()
  })
  
  observeEvent(input$loadButton_example, {
    print("example_button")
    data$master_tbl <- read.table(example_masterfile, header=T, sep="\t", quote="", colClasses=c(Bin="character"))
    data$contig_tbl <- read.table(example_contigfile, header=T, sep="\t", quote="", colClasses=c(Bin="character"))
    data$masterfile <- paste(example_masterfile, '(example)')
    data$contigfile <- paste(example_contigfile, '(example)')
    updateTables()
  })


  
  
  
  ######################  Files summary ############################################
  
  
  output$filesSummary <- renderUI({
    print("files summary")
    
    
    #Calculate some numbers
    metaG_count <- 0
    metaP_count <- 0
    metaT_count <- 0
    cazyme_count <- 0
    kegg_count <- 0
    mag_count <- 0
    if (!is.null(data$masterfile)) {
      metaG_count <- nrow(data$master_tbl)
      if (do_metaP)
        metaP_count <- nrow(data$master_tbl %>% filter(metaP_expressed == TRUE))
      if (do_metaT)
        metaT_count <- nrow(data$master_tbl %>% filter(metaT_expressed == TRUE))
      cazyme_count <- nrow(data$cazymes)
      kegg_count <- nrow(data$kegg_proteins)
      mag_count <- nrow(data$contig_tbl %>% select(Bin) %>% distinct())
    }
    
    
    HTML(paste(
      paste("<b>Masterfile:</b><code>", data$masterfile,"</code>"),
      paste("<b>Contigfile:</b><code>", data$contigfile,"</code>"),
      paste("<b>MAGs:</b><code>", mag_count,"</code>"),
      paste("<b>Genes:</b><code>", metaG_count,"</code>"),
      paste("<b>Transcripts quantified:</b><code>", metaT_count,"</code>"),
      paste("<b>Proteins quantified:</b><code>", metaP_count,"</code>"),
      paste("<b>CAZymes:</b><code>", cazyme_count,"</code>"),
      paste("<b>KEGG-annotated genes:</b><code>", kegg_count,"</code>"),
      paste("<br/>", ifelse(mag_count > 49 | metaG_count > 100000, "<h3><b>NB: You may experience slowness in the App with this many MAGs/genes. It will work, just slow.</b></h3>", "")),
      sep = '<br/>'))
  })
  

    

  ######################  CAZy Plots ###############################################
  
  #Title
  output$cazyCaption <- renderText({ input$cazyGroup })
  
  #Legend
  output$cazyPlot_legend <- renderPlot({

    group <- input$cazyGroup
    dataSubset <- data$cazymes_expanded[data$cazymes_expanded$CAZyText == group,]
    
    colors_to_use <- global_swatch[unique(dataSubset$Bin)]

    tmp_plot <- ggplot(dataSubset) + geom_bar(mapping = aes(x=CAZyClassNumber, fill = Bin)) + scale_fill_manual(values=colors_to_use) + theme(legend.position='bottom')
    legend <- get_legend(tmp_plot)

    grid.newpage()
    grid.draw(legend)
  }, height = function() {
    round(length(data$cazymes_expanded %>%
                   filter(CAZyText == input$cazyGroup) %>%
                   select(Bin) %>% distinct() %>% deframe()) * 2.3 + 30)
  })
  
  #MetaG
  output$cazyPlot_MG <- renderPlot({
    
    showNotification("Loading...", type = "message", duration = 2)

    group <- input$cazyGroup
    dataSubset <- data$cazymes_expanded[data$cazymes_expanded$CAZyText == group,]
    dataSubset <- dataSubset %>% mutate(CAZyClassNumber = factor(CAZyClassNumber, levels = c(1:300, 'cohesin', 'dockerin', 'SLH')))

    plot <- ggplot(dataSubset) +
       geom_bar(mapping = aes(x=CAZyClassNumber, fill = Bin)) +
       scale_fill_manual(values=global_swatch) +
       scale_y_continuous(breaks = function(x) pretty(x, n=3)) +
       labs(title = "", y="Gene count", x="CAZy class number") +
       theme(legend.position='none', axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
    
    ggsave("cazyPlot_MG.pdf", plot)
    plot
    
  })
  output$download_cazyPlot_MG  <- downloadHandler(
    filename = function() {"cazyPlot_MG.pdf"},
    content = function(file) {file.copy("cazyPlot_MG.pdf", file, overwrite=TRUE)}
  )
  
  #MetaT
  output$cazyPlot_MT <- renderPlot({
    print("CAZyPlot_MT rendering")
    group <- input$cazyGroup
    
    if (!do_metaT) {
      print("returning...")
      return()
    }
    
    dataSubset <- data$cazymes_expanded %>%
      select(-contains("metaP_LFQ")) %>%
      filter(CAZyText == group) %>%
      pivot_longer(cols = contains("metaT_TPM"), names_to="Sample", values_to="Value") %>%
      filter(Sample == input$sample_MT)

    dataSubset <- dataSubset %>%
      group_by(CAZyClassNumber, Bin) %>%
      summarise(y_axis = sum(Value, na.rm = TRUE))

    dataSubset <- dataSubset %>% mutate(CAZyClassNumber = factor(CAZyClassNumber, levels = c(1:300, 'cohesin', 'dockerin', 'SLH')))

    plot <- ggplot(dataSubset) +
      geom_col(mapping = aes(x=CAZyClassNumber, fill = Bin, y=y_axis)) +
      scale_fill_manual(values=global_swatch) +
      scale_y_continuous(breaks = function(x) pretty(x, n=3)) +
      labs(title = "", y="Transcripts per million (tpm)", x="CAZy class number") +
      theme(legend.position='none', axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
    
    ggsave("cazyPlot_MT.pdf", plot)
    plot
  })
  output$download_cazyPlot_MT  <- downloadHandler(
    filename = function() {"cazyPlot_MT.pdf"},
    content = function(file) {file.copy("cazyPlot_MT.pdf", file, overwrite=TRUE)}
  )

  #MetaP
  output$cazyPlot_MP <- renderPlot({
    print("CAZyPlot_MP rendering")
    group <- input$cazyGroup

    if (!do_metaP) {
      print("returning...")
      return()
    }
    
    
    dataSubset <- data$cazymes_expanded %>%
      select(-contains("metaT_TPM")) %>%
      filter(CAZyText == group) %>%
      pivot_longer(cols = contains("metaP_LFQ"), names_to="Sample", values_to="Value") %>%
      filter(Sample == input$sample_MP)

    if (input$yaxis_MP == "Summed LFQ") {
      dataSubset <- dataSubset %>%
        group_by(CAZyClassNumber, Bin) %>%
        summarise(y_axis = sum(Value, na.rm = TRUE))
    }
    if (input$yaxis_MP == "Protein count") {
      dataSubset <- dataSubset %>%
        group_by(CAZyClassNumber, Bin) %>%
        summarise(y_axis = sum(!is.na(Value)))
    }

    dataSubset <- dataSubset %>% mutate(CAZyClassNumber = factor(CAZyClassNumber, levels = c(1:300, 'cohesin', 'dockerin', 'SLH')))


    plot <- ggplot(dataSubset) +
      geom_col(mapping = aes(x=CAZyClassNumber, fill = Bin, y=y_axis)) +
      scale_fill_manual(values=global_swatch) +
      scale_y_continuous(breaks = function(x) pretty(x, n=3)) +
      labs(title = "", y=input$yaxis_MP, x="CAZy class number") +
      theme(legend.position='none', axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
    
    ggsave("cazyPlot_MP.pdf", plot)
    plot
    
  })
  output$download_cazyPlot_MP  <- downloadHandler(
    filename = function() {"cazyPlot_MP.pdf"},
    content = function(file) {file.copy("cazyPlot_MP.pdf", file, overwrite=TRUE)}
  )

  #CAZy-heatmaps, updates based on selected bin
  output$cazyPlot_heatmap_MT <- renderPlot({
    print("CAZyPlot_heatmap_MT rendering")
    
    group <- input$cazyGroup

    if (!do_metaT) {
      print("returning...")
      return()
    }
    
    
    selected_data <- data$cazymes_expanded %>%
      select(-contains("metaP_LFQ")) %>%
      filter(CAZyText == group, Bin == input$cazyPlot_heatmap_binSelector) %>%
      pivot_longer(cols = contains("metaT_TPM"), names_to = "Sample", values_to = "Value") %>%
      group_by(CAZyClassNumber, Sample) %>%
      summarise(Value = sum(Value, na.rm = TRUE))
    
    n_cc <- selected_data %>% select(CAZyClassNumber) %>% distinct() %>% deframe()

    if (length(n_cc) == 0) {
      ggplot() + annotate("text", x=4, y=25, size=8, label="Not enough CAZy data for the selected bin") + theme_void()
    } else {
      
      #Sort numbers alphabetically
      selected_data <- selected_data %>% mutate(CAZyClassNumber = factor(CAZyClassNumber, levels = c(1:300, 'cohesin', 'dockerin', 'SLH')))

      plot <- ggplot(selected_data) +
        geom_tile(aes(y=CAZyClassNumber, x=Sample, fill=Value), color="white") +
        labs(title = paste0("Abundance of ", group,", summed transcripts per million (tpm)")) +
        theme(legend.position = 'top', axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
        guides(fill = guide_colourbar(title="", barwidth = 20, barheight = 0.5)) +
        scale_y_discrete(position = "right")

      ggsave("cazyPlot_MT_heatmap.pdf", plot)
      plot
      
    }
  },width = function() {
    ncol(data$cazymes_expanded %>% select(contains("metaT_TPM"))) * 50 + 200
    
  }, height = function() {
    round(length(data$cazymes_expanded %>%
            filter(CAZyText == input$cazyGroup, Bin == input$cazyPlot_heatmap_binSelector) %>%
            select(CAZyClassNumber) %>% distinct() %>% deframe()) * 15 + 250)
  })
  output$download_cazyPlot_MT_heatmap  <- downloadHandler(
    filename = function() {"cazyPlot_MT_heatmap.pdf"},
    content = function(file) {file.copy("cazyPlot_MT_heatmap.pdf", file, overwrite=TRUE)}
  )

  
  output$cazyPlot_heatmap_MP <- renderPlot({
    print("CAZyPlot_heatmap_MP rendering")
    
    group <- input$cazyGroup
    
    if (!do_metaP) {
      print("returning...")
      return()
    }
    
    
    selected_data <- data$cazymes_expanded %>%
      select(-contains("metaT_TPM")) %>%
      filter(CAZyText == group, Bin == input$cazyPlot_heatmap_binSelector) %>%
      pivot_longer(cols = contains("metaP_LFQ"), names_to = "Sample", values_to = "Value") %>%
      group_by(CAZyClassNumber, Sample) %>%
      summarise(Value = sum(Value, na.rm = TRUE))
    
    n_cc <- selected_data %>% select(CAZyClassNumber) %>% distinct() %>% deframe()

    if (length(n_cc) == 0) {
      ggplot() + annotate("text", x=4, y=25, size=8, label="Not enough CAZy data for the selected bin") + theme_void()
    } else {
      
      #Sort numbers alphabetically
      selected_data <- selected_data %>% mutate(CAZyClassNumber = factor(CAZyClassNumber, levels = c(1:300, 'cohesin', 'dockerin', 'SLH')))

      plot <- ggplot(selected_data) +
        geom_tile(aes(y=CAZyClassNumber, x=Sample, fill=Value), color="white") +
        labs(title = paste0("Abundance of ", group,", summed protein LFQ")) +
        theme(legend.position = 'top', axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
        guides(fill = guide_colourbar(title="", barwidth = 20, barheight = 0.5)) +
        scale_y_discrete(position = "right")
      
      ggsave("cazyPlot_MP_heatmap.pdf", plot)
      plot
    }
  },width = function() {
    ncol(data$cazymes_expanded %>% select(contains("metaP_LFQ"))) * 50 + 200
    
  },height = function() {
    round(length(data$cazymes_expanded %>%
                   filter(CAZyText == input$cazyGroup, Bin == input$cazyPlot_heatmap_binSelector) %>%
                   select(CAZyClassNumber) %>% distinct() %>% deframe()) * 15 + 250)
  })
  output$download_cazyPlot_MP_heatmap  <- downloadHandler(
    filename = function() {"cazyPlot_MP_heatmap.pdf"},
    content = function(file) {file.copy("cazyPlot_MP_heatmap.pdf", file, overwrite=TRUE)}
  )


  
  
  
  

  
  
  
  
  ######################  KEGG Plots ###############################################
  
  output$keggCaption <- renderText({ input$keggGroup })

  #Plots, updates on selected function in combo-box
  output$keggPlot_legend <- renderPlot({

    colors_to_use <- global_swatch[unique(data$kegg_proteins_expanded$Bin)]
    
    tmp_plot <- ggplot(data$kegg_proteins_expanded) + geom_bar(mapping = aes(x=KEGG, fill = Bin)) + scale_fill_manual(values=colors_to_use) + theme(legend.position='bottom')
    legend <- get_legend(tmp_plot)

    grid.newpage()
    grid.draw(legend)
  }, height = function() {
    round(length(data$kegg_proteins_expanded %>%
                   filter(Parent.term == input$keggGroup) %>%
                   select(Bin) %>% distinct() %>% deframe()) * 2.3 + 60)
  })

  #KEGGPlot - MG
  output$keggPlot_MG <- renderPlot({
    
    showNotification("Loading...", type = "message", duration = 2)
    

    dataSubset <- data$kegg_proteins_expanded %>% filter(Parent.term == input$keggGroup)

    plot <- ggplot(dataSubset) +
      geom_bar(mapping = aes(x=KEGG, fill = Bin)) +
      scale_fill_manual(values=global_swatch) +
      scale_y_continuous(breaks = function(x) pretty(x, n=3)) +
      labs(title = "", y="Gene count", x="KEGG Pathway") +
      theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    ggsave("keggPlot_MG.pdf", plot)
    plot
  })
  output$download_keggPlot_MG  <- downloadHandler(
    filename = function() {"keggPlot_MG.pdf"},
    content = function(file) {file.copy("keggPlot_MG.pdf", file, overwrite=TRUE)}
  )

  #KEGGPlot - MT
  output$keggPlot_MT <- renderPlot({
    print("KeggPlot_MT rendering")
    
    group <- input$keggGroup
    
    if (!do_metaT) {
      print("returning...")
      return()
    }
    
    dataSubset <- data$kegg_proteins_expanded %>%
      select(-contains("metaP_LFQ")) %>%
      filter(Parent.term == input$keggGroup) %>%
      pivot_longer(cols = contains("metaT_TPM"), names_to="Sample", values_to="Value") %>%
      filter(Sample == input$kegg_sample_MT)

    dataSubset <- dataSubset %>%
      group_by(KEGG, Bin) %>%
      summarise(y_axis = sum(Value, na.rm = TRUE))

    plot <- ggplot(dataSubset) +
      geom_col(mapping = aes(x=KEGG, fill = Bin, y=y_axis)) +
      scale_fill_manual(values=global_swatch) +
      scale_y_continuous(breaks = function(x) pretty(x, n=3)) +
      labs(title = "", y="Transcripts per million (tpm)", x="KEGG Pathway") +
      theme(legend.position='none', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    ggsave("keggPlot_MT.pdf", plot)
    plot
  })
  output$download_keggPlot_MT  <- downloadHandler(
    filename = function() {"keggPlot_MT.pdf"},
    content = function(file) {file.copy("keggPlot_MT.pdf", file, overwrite=TRUE)}
  )


  #KEGGPlot - MP
  output$keggPlot_MP <- renderPlot({
    print("KeggPlot_MP rendering")
    
    group <- input$keggGroup
    
    if (!do_metaP) {
      print("returning...")
      return()
    }
    
    dataSubset <- data$kegg_proteins_expanded %>%
      select(-contains("metaT_TPM")) %>%
      filter(Parent.term == input$keggGroup) %>%
      pivot_longer(cols = contains("metaP_LFQ"), names_to="Sample", values_to="Value") %>%
      filter(Sample == input$kegg_sample_MP)

    if (input$kegg_yaxis_MP == "Summed LFQ") {
      dataSubset <- dataSubset %>%
        group_by(KEGG, Bin) %>%
        summarise(y_axis = sum(Value, na.rm = TRUE))
    }
    if (input$kegg_yaxis_MP == "Protein count") {
      dataSubset <- dataSubset %>%
        group_by(KEGG, Bin) %>%
        summarise(y_axis = sum(!is.na(Value)))
    }

    plot <- ggplot(dataSubset) +
      geom_col(mapping = aes(x=KEGG, fill = Bin, y=y_axis)) +
      scale_fill_manual(values=global_swatch) +
      scale_y_continuous(breaks = function(x) pretty(x, n=3)) +
      labs(title = "", y=input$kegg_yaxis_MP, x="KEGG Pathway") +
      theme(legend.position='none', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    ggsave("keggPlot_MP.pdf", plot)
    plot
  })
  output$download_keggPlot_MP  <- downloadHandler(
    filename = function() {"keggPlot_MP.pdf"},
    content = function(file) {file.copy("keggPlot_MP.pdf", file, overwrite=TRUE)}
  )


  #Need a proxy to manipulate the DT table
  #proxy_keggTable <- dataTableProxy('keggTable')

  #Data subset generated by selecting a parent function
  keggTable_dataSubset <- reactive({
    data$kegg_proteins_expanded %>%
      filter(Parent.term == input$keggGroup) %>%
      select(KEGG, Name) %>%
      distinct()
    
    #Now the parent term changed. Remove any selection from pathway-table
    #proxy_keggTable %>% selectRows(NULL)
  })

  #Table, updates on selected function, interactive by click
  output$keggTable <- DT::renderDataTable(
    keggTable_dataSubset(),
    server = FALSE,
    selection = list(mode = 'single', target = "row", selected = NULL),
    class = 'compact hover',
    extensions = 'Buttons',
    options = list(
      dom = "Blfrtip",
      buttons = c('copy', 'csv')),
    rownames = FALSE,
    caption = paste0('This table lists the KEGG pathways within the domain "', input$keggGroup, '" and found in the dataset')
    )




  #KEGG-heatmaps, updates based on selected row in table, and on selected bin
  output$keggPlot_heatmap_MT <- renderPlot({
    print("KeggPlot_heatmap_MT rendering")
    
    s = input$keggTable_rows_selected
    
    if (!do_metaT) {
      print("returning...")
      return()
    }
    
    if (length(s)) {
      selected_kegg <- as.character(keggTable_dataSubset()[s,] %>% select(KEGG))

      # selected_data <- data$kegg_proteins_expanded %>%
      #   select(-contains("metaP_LFQ")) %>%
      #   filter(Parent.term == input$keggGroup, KEGG == selected_kegg, Bin == input$keggPlot_heatmap_binSelector, EC != "") %>%
      #   pivot_longer(cols = contains("metaT_TPM"), names_to = "Sample", values_to = "Value") %>%
      #   group_by(EC, Sample) %>%
      #   summarise(Value = sum(Value, na.rm = TRUE))
      selected_data <- data$kegg_proteins_expanded %>%
        select(-contains("metaP_LFQ")) %>%
        filter(Parent.term == input$keggGroup, KEGG == selected_kegg, Bin == input$keggPlot_heatmap_binSelector, KO != "") %>%
        mutate(Enzyme = ifelse(EC=="", KO, paste0(KO, " (", EC, ")"))) %>%
        pivot_longer(cols = contains("metaT_TPM"), names_to = "Sample", values_to = "Value") %>%
        group_by(Enzyme, Sample) %>%
        summarise(Value = sum(Value, na.rm = TRUE))
      
      n_ec <- selected_data %>% select(Enzyme) %>% distinct() %>% deframe()

      if (length(n_ec) == 0) {
        ggplot() + annotate("text", x=4, y=25, size=8, label="Not enough KEGG data for the selected bin") + theme_void()
      } else {

        plot <- ggplot(selected_data) +
          geom_tile(aes(y=Enzyme, x=Sample, fill=Value), color="white") +
          labs(title = paste0("Abundance of enzymes in pathway ", selected_kegg,", Transcripts per million (tpm)")) +
          theme(legend.position = 'top', axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
          guides(fill = guide_colourbar(title="", barwidth = 20, barheight = 0.5)) +
          scale_y_discrete(position = "right")
        
        ggsave("keggPlot_MT_heatmap.pdf", plot)
        plot
        
      }

    }
  },width = function() {
    ncol(data$kegg_proteins_expanded %>% select(contains("metaT_TPM"))) * 30 + 550
    
  },height = function() {
    round(length(data$kegg_proteins_expanded %>%
                   filter(Parent.term == input$keggGroup, KEGG == as.character(keggTable_dataSubset()[input$keggTable_rows_selected,] %>% select(KEGG)), Bin == input$keggPlot_heatmap_binSelector) %>%
                   select(KO) %>% distinct() %>% deframe()) * 15 + 250)
  })
  output$download_keggPlot_MT_heatmap  <- downloadHandler(
    filename = function() {"keggPlot_MT_heatmap.pdf"},
    content = function(file) {file.copy("keggPlot_MT_heatmap.pdf", file, overwrite=TRUE)}
  )

  
  
  output$keggPlot_heatmap_MP <- renderPlot({
    print("KeggPlot_heatmap_MP rendering")
    
    s = input$keggTable_rows_selected
    
    if (!do_metaP) {
      print("returning...")
      return()
    }
    
    if (length(s)) {
      selected_kegg <- as.character(keggTable_dataSubset()[s,] %>% select(KEGG))

      # selected_data <- data$kegg_proteins_expanded %>%
      #   select(-contains("metaT_TPM")) %>%
      #   filter(Parent.term == input$keggGroup, KEGG == selected_kegg, Bin == input$keggPlot_heatmap_binSelector, EC != "") %>%
      #   pivot_longer(cols = contains("metaP_LFQ"), names_to = "Sample", values_to = "Value") %>%
      #   group_by(EC, Sample) %>%
      #   summarise(Value = sum(Value, na.rm = TRUE))
      selected_data <- data$kegg_proteins_expanded %>%
        select(-contains("metaT_TPM")) %>%
        filter(Parent.term == input$keggGroup, KEGG == selected_kegg, Bin == input$keggPlot_heatmap_binSelector, KO != "") %>%
        mutate(Enzyme = ifelse(EC=="", KO, paste0(KO, " (", EC, ")"))) %>%
        pivot_longer(cols = contains("metaP_LFQ"), names_to = "Sample", values_to = "Value") %>%
        group_by(Enzyme, Sample) %>%
        summarise(Value = sum(Value, na.rm = TRUE))
      
      n_ec <- selected_data %>% select(Enzyme) %>% distinct() %>% deframe()

      if (length(n_ec) == 0) {
        ggplot() + annotate("text", x=4, y=25, size=8, label="Not enough KEGG data for the selected bin") + theme_void()
      } else {
        
        plot <- ggplot(selected_data) +
          geom_tile(aes(y=Enzyme, x=Sample, fill=Value), color="white") +
          labs(title = paste0("Abundance of enzymes in pathway ", selected_kegg,", Protein LFQ")) +
          theme(legend.position = 'top', axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
          guides(fill = guide_colourbar(title="", barwidth = 20, barheight = 0.5)) +
          scale_y_discrete(position = "right")
        
        ggsave("keggPlot_MP_heatmap.pdf", plot)
        plot
        
      }
    }
  },width = function() {
    ncol(data$kegg_proteins_expanded %>% select(contains("metaP_LFQ"))) * 30 + 550
    
  },height = function() {
    round(length(data$kegg_proteins_expanded %>%
                   filter(Parent.term == input$keggGroup, KEGG == as.character(keggTable_dataSubset()[input$keggTable_rows_selected,] %>% select(KEGG)), Bin == input$keggPlot_heatmap_binSelector) %>%
                   select(KO) %>% distinct() %>% deframe()) * 15 + 250)
  })
  output$download_keggPlot_MP_heatmap  <- downloadHandler(
    filename = function() {"keggPlot_MP_heatmap.pdf"},
    content = function(file) {file.copy("keggPlot_MP_heatmap.pdf", file, overwrite=TRUE)}
  )
    

  #Colored KEGG pathway
  output$keggPathway <- renderImage({
    s = input$keggTable_rows_selected
    if (length(s) && input$keggplot_check_download_kegg) {

      selected_kegg <- as.character(keggTable_dataSubset()[s,] %>% select(KEGG))

      #Different coloring of the KEGG pathway depending on data set chosen (MG, MT, MP)

      if (input$keggplot_annotation_radio == "MetaG") {
        
        #Since 230314, I moved from EC to KO as these are more stable and "mappable" at KEGG
        #Hence better results and less errors with KEGG REST.
        
        # ecs <- data$kegg_proteins_expanded %>%
        #   filter(Parent.term == input$keggGroup, KEGG == selected_kegg, Bin == input$keggPlot_heatmap_binSelector, EC != "") %>%
        #   select(EC) %>%
        #   separate_rows(EC, sep = ";") %>%
        #   distinct() %>%
        #   deframe()

        kos <- data$kegg_proteins_expanded %>%
          filter(Parent.term == input$keggGroup, KEGG == selected_kegg, Bin == input$keggPlot_heatmap_binSelector, KO != "") %>%
          select(KO) %>%
          separate_rows(KO, sep = ";") %>%
          distinct() %>%
          filter(KO %in% k_numbers_connected_to_pathways) %>%
          deframe()
        

        if (length(kos) == 0)
          return(list(src = tempfile(fileext='.png'), alt= "Not enough KEGG data for the selected bin"))
        
        #pathway <- gsub('ko', 'ec', selected_kegg, fixed = TRUE) #Needed when using ECs. The ko-maps are only for K-numbers

        outfile <- tempfile(fileext='.png')
        image_url <- ''
        
        tryCatch(
          expr = {
            pathway_id <- paste0("path:", selected_kegg)
            print(pathway_id)
            print(kos)
            
            image_url <- mark.pathway.by.objects(pathway_id, kos)
            
            print(image_url)
            
            download.file(image_url, outfile, method = "auto", mode = "wb")
            output$ERRORtext_KEGGmap <- renderText({ "" })
          },
          error = function(e) {
            print("An error occured in output$keggPathway while contacting KEGG")
            print(e)
            
            t <- paste("<h4 style='color:red;'>An error occured while contacting KEGG. Perhaps the KEGG REST library needs an update. Please contact the developers.<br><br>Here is the error:<br>", e)
            output$ERRORtext_KEGGmap <- renderText({ t })
          },
          warning = function(e) {
            print("A warning occured in output$keggPathway while contacting KEGG")
            print(e)
            t <- paste("<h4 style='color:red;'>A warning occured while contacting KEGG. Perhaps the KEGG REST library needs an update. Please contact the developers.<br><br>Here is the warning:<br>", e)
            output$ERRORtext_KEGGmap <- renderText({ t })
          }
        )

        #Return a list with the image.
        return(list(src = outfile, alt= image_url, style="max-width: 100%;"))
      }
      else {
        #metaT or metaP
        
        selected_data <- NULL
        
        
        #Since 230314, I moved from EC to KO as these are more stable and "mappable" at KEGG
        #Hence better results and less errors with KEGG REST.
        
        
        if (do_metaT && input$keggplot_annotation_radio == "MetaT") {
          # selected_data <- data$kegg_proteins_expanded %>%
          #   select(-contains("metaP_LFQ")) %>%
          #   filter(Parent.term == input$keggGroup, KEGG == selected_kegg, Bin == input$keggPlot_heatmap_binSelector) %>%
          #   pivot_longer(cols = contains("metaT_TPM"), names_to = "Sample", values_to = "Value") %>%
          #   separate_rows(EC, sep = ";") %>%
          #   group_by(EC) %>%
          #   summarise(Average = mean(Value, na.rm = TRUE)) %>%
          #   mutate(Average = replace_na(Average, 0))

          selected_data <- data$kegg_proteins_expanded %>%
            select(-contains("metaP_LFQ")) %>%
            filter(Parent.term == input$keggGroup, KEGG == selected_kegg, Bin == input$keggPlot_heatmap_binSelector) %>%
            pivot_longer(cols = contains("metaT_TPM"), names_to = "Sample", values_to = "Value") %>%
            separate_rows(KO, sep = ";") %>%
            filter(KO %in% k_numbers_connected_to_pathways) %>%
            group_by(KO) %>%
            summarise(Average = mean(Value, na.rm = TRUE)) %>%
            mutate(Average = replace_na(Average, 0))
        }
        if (do_metaP && input$keggplot_annotation_radio == "MetaP") {
          # selected_data <- data$kegg_proteins_expanded %>%
          #   select(-contains("metaT_TPM")) %>%
          #   filter(Parent.term == input$keggGroup, KEGG == selected_kegg, Bin == input$keggPlot_heatmap_binSelector) %>%
          #   pivot_longer(cols = contains("metaP_LFQ"), names_to = "Sample", values_to = "Value") %>%
          #   separate_rows(EC, sep = ";") %>%
          #   group_by(EC) %>%
          #   summarise(Average = mean(Value, na.rm = TRUE)) %>%
          #   mutate(Average = replace_na(Average, 0))

          selected_data <- data$kegg_proteins_expanded %>%
            select(-contains("metaT_TPM")) %>%
            filter(Parent.term == input$keggGroup, KEGG == selected_kegg, Bin == input$keggPlot_heatmap_binSelector) %>%
            pivot_longer(cols = contains("metaP_LFQ"), names_to = "Sample", values_to = "Value") %>%
            separate_rows(KO, sep = ";") %>%
            filter(KO %in% k_numbers_connected_to_pathways) %>%
            group_by(KO) %>%
            summarise(Average = mean(Value, na.rm = TRUE)) %>%
            mutate(Average = replace_na(Average, 0))
        }
        
        if (!is.null(selected_data)) {
        
          #ecs <- selected_data %>% select(EC) %>% deframe()
          kos <- selected_data %>% select(KO) %>% deframe()
          numbers <- selected_data %>% select(Average) %>% deframe()

          if (length(kos) == 0)
            return(list(src = tempfile(fileext='.png'), alt= "Not enough KEGG data for the selected bin"))
          

          #Make a palette with 7 colors and their corresponding intensities
          low <- '#fff2a9'
          high <- '#a7372b'
          col <- colorRampPalette(c(low, high)) (7) #7 steps are OK for KEGG visualization
          lowest_number <- min(numbers)
          highest_number <- max(numbers)
          added = FALSE
          if (highest_number == lowest_number) {
            highest_number = highest_number + 50 #to get a gradient
            numbers <- c(numbers, highest_number) #need also this high number to get the map2color correct, will remove before calling KEGG-Map
            added = TRUE
          }
          col_numbers <- round(seq.int(lowest_number, highest_number, length.out = 7))
          
          #Make legend color palette image
          png(filename = "KEGG_legend.png", width = 300, height = 100)
          par(mar=c(4,1,1,1))
          image(col_numbers, 1, as.matrix(1:length(col)), xlab="", ylab="", yaxt="n")
          dev.off()
          
          #Map all intesities to their respective colors and download the map
          numbers <- map2color(numbers, col)
  
          #pathway <- gsub('ko', 'ec', selected_kegg, fixed = TRUE) #Needed when using ECs. The ko-maps are only for K-numbers
          pathway <- selected_kegg
          
          if (added) numbers <- numbers[-length(numbers)] #remove again

          image_url <- ''
          success <- TRUE
          
          #Contact KEGG to retrieve image
          tryCatch(
            expr = {
              image_url <- color.pathway.by.objects(paste0("path:", pathway), kos, rep_len('black', length.out = length(kos)), numbers)
              download.file(image_url, "KEGG_map.png", method = "auto", mode = "wb")
              
              
              #Make a composite image with the two
              composite <- image_append(c(image_read("KEGG_map.png"), image_read("KEGG_legend.png")), stack = TRUE)
              image_write(composite, path = paste0(wd, "/KEGG_map_legend.png"), format = "png")
              
              #No error
              output$ERRORtext_KEGGmap <- renderText({ "" })
            },
            error = function(e) {
              print("An error occured in output$keggPathway while contacting KEGG")
              print(e)
              
              t <- paste("<h4 style='color:red;'>An error occured while contacting KEGG. Perhaps the KEGG REST library needs an update. Please contact the developers.<br><br>Here is the error:<br>", e)
              output$ERRORtext_KEGGmap <- renderText({ t })
              
              success <- FALSE
            },
            warning = function(e) {
              print("A warning occured in output$keggPathway while contacting KEGG")
              print(e)
              t <- paste("<h4 style='color:red;'>A warning occured while contacting KEGG. Perhaps the KEGG REST library needs an update. Please contact the developers.<br><br>Here is the warning:<br>", e)
              output$ERRORtext_KEGGmap <- renderText({ t })
              
              success <- FALSE
            }
          )
          
          if (success) {
            outfile <- paste0(wd, "/KEGG_map_legend.png")
          } else {
            outfile <- tempfile(fileext='.png')
          }
          
          #Return a list with the image.
          return(list(src = outfile, alt= image_url, style="max-width: 100%;"))
        }
      }


    } else {

      outfile <- tempfile(fileext='.png')
      return(list(src = outfile, alt= "No KEGG file available"))
    }
  }, deleteFile = TRUE)




  
  
  
  ######################  MAG Plots ###############################################
  
  
  #The data subset needs to be in a reactive environment as data$contig_tbl is reactive
  MAGTable_dataSubset <- reactive({

    datasubset <- data$contig_tbl %>%
                group_by(Bin) %>%
                summarise(`Contig count` = sum(!is.na(Contig)),
                          `Avg. Coverage` = mean(Coverage, na.rm = TRUE)) %>%
                left_join(data$master_tbl %>%
                            select(c('Accn', 'Bin')) %>%
                            group_by(Bin) %>%
                            summarise(`Gene count` = sum(!is.na(Accn))),
                            by = c("Bin" = "Bin"))
    
    if (do_lineage) {
      datasubset <- datasubset %>%
        left_join(data$contig_tbl %>%
                    select(Bin, Lineage) %>%
                    distinct(),
                  by = c("Bin" = "Bin"))
    }

    if (do_checkM) {
      datasubset <- datasubset %>%
        left_join(data$contig_tbl %>%
                    select(Bin, Completeness, Contamination, Strain.heterogeneity) %>%
                    distinct(),
                  by = c("Bin" = "Bin"))
    }
    
    
    #Filter to the sliders
    datasubset <- datasubset %>% 
      filter(Completeness >= input$MAG_slider_completeness, Contamination <= input$MAG_slider_contamination)

    datasubset
  })

  
  #Static, but interactive table
  output$MAGTable <- DT::renderDataTable(MAGTable_dataSubset(),
                                         server = FALSE,
                                         selection = 'single',
                                         class = 'compact hover',
                                         extensions = 'Buttons',
                                         options = list(
                                           dom = "Blfrtip",
                                           buttons = c('copy', 'csv')),
                                         rownames = FALSE)


  #MAG-plot info string
  output$MAGplot_info <- renderText({
    
    bin_count <- length(unique(data$contig_tbl$Bin))
    contig_count <- length(unique(data$contig_tbl$Contig))

    paste0("There are ",bin_count ," MAGs and ", contig_count, " contigs in your data")
  })
  
  
  
  #Plot, with highlighting from table
  output$MAGPlot <- renderPlot({
    print("rendering MAG plot")
    
    showNotification("Loading...", type = "message", duration = 2)
    

    #Plot shows the MAGs from the table, i.e. affected by the filtering options
    bins_to_show <- MAGTable_dataSubset() %>% select(Bin) %>% distinct() %>% deframe()
    
    dataSubset <- data$contig_tbl %>%
      rename(`Contig length` = Length) %>%
      filter(Bin %in% bins_to_show)
    
    moreThan50 <- length(unique(dataSubset$Bin)) > 49
    
    #Limit the plot to the 50 most abundant MAGs
    if (moreThan50) {
      mostAbundantMAGS <- MAGTable_dataSubset() %>% slice_max(`Avg. Coverage`, n = 50) %>% select(Bin) %>%  distinct() %>% deframe()
      dataSubset <- dataSubset %>% filter(Bin %in% mostAbundantMAGS)
      
      print("Too many MAGs for plotting. Showing only top50.")
    }
    
    colors_to_use <- global_swatch[unique(dataSubset$Bin)]
    
    plot <- ggplot(dataSubset, aes(x=GC, y=Coverage, size=`Contig length`, color=Bin)) +
      geom_point(alpha=0.7) +
      scale_color_manual(values=colors_to_use) +
      labs(title="", x="%GC", y="Coverage") +
      scale_y_log10()

    s = input$MAGTable_rows_selected

    if (length(s)) {
      selected_bin <- as.character(MAGTable_dataSubset()[s,] %>% select(Bin))
      selected_data <- dataSubset %>%
        filter(Bin %in% c(selected_bin))

      #Could be that the selected bin is not drawn if the plot is limited to only the most abundant
      if (length(selected_data))
        plot <- plot + geom_point(data = selected_data, mapping = aes(x=GC, y=Coverage, size=`Contig length`, color=Bin), alpha=0.9, shape=21, fill="white", stroke=2)
    }

    ggsave("MAG_GC_coverage.pdf", plot)
    plot
  })

  output$downloadMAGPlot <- downloadHandler(
    filename = function() {"MAG_GC_coverage.pdf"},
    content = function(file) {file.copy("MAG_GC_coverage.pdf", file, overwrite=TRUE)}
  )


  

  #Heatmap - MAG abundance MT
  output$MAGabundancePlot_heatmap_MT <- renderPlot({
    
    selected_data <- data$master_tbl %>%
      select(Accn, Bin, contains("metaT_TPM")) %>%
      pivot_longer(cols = contains("metaT_TPM"), names_to = "Sample", values_to = "Value") %>%
      group_by(Bin, Sample) %>%
      summarise(MAG_abundance = sum(Value, na.rm = TRUE))
    
    if (!do_metaT) return()

    n_cc <- selected_data %>% select(Bin) %>% distinct() %>% deframe()
    
    if (length(n_cc) == 0) {
      ggplot() + annotate("text", x=4, y=25, size=8, label="Not enough MAG data for the selected bin") + theme_void()
    } else {
      
      #Sort numbers alphabetically
      #selected_data <- selected_data %>% mutate(CAZyClassNumber = factor(CAZyClassNumber, levels = c(1:300, 'cohesin', 'dockerin', 'SLH')))
      
      plot <- ggplot(selected_data) +
        geom_tile(aes(y=Bin, x=Sample, fill=MAG_abundance), color="white") +
        labs(title = "MAG Abundance per sample, summed transcripts per million (tpm)") +
        theme(legend.position = 'top', axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
        guides(fill = guide_colourbar(title="", barwidth = 20, barheight = 0.5)) +
        scale_y_discrete(position = "right")
      
      ggsave("MAG_abundancePlot_MT_heatmap.pdf", plot)
      plot
    }
  },width = function() {
    ncol(data$master_tbl %>% select(contains("metaT_TPM"))) * 50 + 200
    
  },height = function() {
    round(length(data$master_tbl %>%
                   select(Bin) %>% distinct() %>% deframe()) * 15 + 250)
  })
  output$download_MAGabundancePlot_MT_heatmap  <- downloadHandler(
    filename = function() {"MAG_abundancePlot_MT_heatmap.pdf"},
    content = function(file) {file.copy("MAG_abundancePlot_MT_heatmap.pdf", file, overwrite=TRUE)}
  )
    

  #Heatmap - MAG abundance MP
  output$MAGabundancePlot_heatmap_MP <- renderPlot({
    
    selected_data <- data$master_tbl %>%
      select(Accn, Bin, contains("metaP_LFQ")) %>%
      pivot_longer(cols = contains("metaP_LFQ"), names_to = "Sample", values_to = "Value") %>%
      group_by(Bin, Sample) %>%
      summarise(MAG_abundance = sum(Value, na.rm = TRUE))
    
    if (!do_metaP) return()
    
    n_cc <- selected_data %>% select(Bin) %>% distinct() %>% deframe()
    
    if (length(n_cc) == 0) {
      ggplot() + annotate("text", x=4, y=25, size=8, label="Not enough MAG data for the selected bin") + theme_void()
    } else {
      
      #Sort numbers alphabetically
      #selected_data <- selected_data %>% mutate(CAZyClassNumber = factor(CAZyClassNumber, levels = c(1:300, 'cohesin', 'dockerin', 'SLH')))
      
      plot <- ggplot(selected_data) +
        geom_tile(aes(y=Bin, x=Sample, fill=MAG_abundance), color="white") +
        labs(title = "MAG Abundance per sample, summed protein LFQ") +
        theme(legend.position = 'top', axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
        guides(fill = guide_colourbar(title="", barwidth = 20, barheight = 0.5)) +
        scale_y_discrete(position = "right")
      
      ggsave("MAG_abundancePlot_MP_heatmap.pdf", plot)
      plot
    }
  },width = function() {
    ncol(data$master_tbl %>% select(contains("metaP_LFQ"))) * 50 + 200
    
  },height = function() {
    round(length(data$master_tbl %>%
                   select(Bin) %>% distinct() %>% deframe()) * 15 + 250)
  })
  output$download_MAGabundancePlot_MP_heatmap  <- downloadHandler(
    filename = function() {"MAG_abundancePlot_MP_heatmap.pdf"},
    content = function(file) {file.copy("MAG_abundancePlot_MP_heatmap.pdf", file, overwrite=TRUE)}
  )
    

  
  
  
  
  ######################  Settings tab #############################################
  
    
  #Theme selector
  output$gg_themeCaption <- renderText({
    print(paste("Settings changed, new ggplot theme:", input$gg_theme))
   
    setThemeAndCalculateColors(input$gg_theme)

    "Remember to redraw graphs for the new theme to take effect!"
    #a("See the the themes here", href="https://github.com/Mikata-Project/ggthemr")
  })
  
  

  
  
  
  
  
    
  ######################  MFC Plots ###############################################
  
  #Data subset generated by selecting a parent function
  mcf_data <- reactiveValues (
    raw_results = NULL
    
    #can I make some pretty-functions within here?
    
  )
  
  
  
  
  observeEvent(input$calculate_mcf_button, {
    print("mcf_button")
    
    #check if running already
    if (mcf_running) {
      showNotification("MCF calculation already in progress", type = "warning")
      return()
    } else {
      showNotification("MCF calculation started...This may take several minutes!", type = "message", duration = NULL)
      mcf_running <<- TRUE
      
      #Do calculation here
      #Can we do asyncronous or will it freeze app?
      
      bin_vs_KO <- data$master_tbl %>%
        select(Bin, KO) %>%
        mutate(KO = na_if(KO, '')) %>%
        filter(!is.na(KO)) %>%
        separate_rows(KO, sep = ";") %>%
        group_by(Bin) %>%
        summarise(K = paste0(KO, collapse = ";"))
      
      
      #Put the results into a reactive object so that the plots will update automatically when this changes
      print("starting mcf...")
      mcf_data$raw_results <- query_genomes_to_modules(as.data.frame(bin_vs_KO), GENOME_ID_COL = "Bin", GENES_COL = "K", META_OUT = T, ADD_OUT = T)
      print("ended mcf...")
      mcf_running <<- FALSE
      
      #MCF-page
      updateSelectInput(session, "mcf_keggGroup", choices = c('No filtering', unique(mcf_data$raw_results$METADATA$CLASS_III)))
      
      
    }
      
  })
  
  output$mcfPlot <- renderPlot({
    print("rendering MCF plot")
    
    if (is.null(mcf_data$raw_results)) return()
    
    
    
    
    mcf_filtered <- as.data.frame(mcf_data$raw_results$MATRIX)
    
    
    #Filter the modules to the selection box of Class_III terms
    if (input$mcf_keggGroup != 'No filtering') {
      mcf_filtered <- as.data.frame(t(mcf_filtered))
      mcf_filtered$CLASS_III <- mcf_data$raw_results$METADATA$CLASS_III
      mcf_filtered <- mcf_filtered %>%
        filter(CLASS_III == input$mcf_keggGroup) %>%
        select(-CLASS_III)
      mcf_filtered <- as.data.frame(t(mcf_filtered))

    }


    #Filter the modules to at least 0.5 mcf in at least one bin
    #EDIT: the cutoff value is now changeable for the user
    cutoff <- input$mcf_module_cutoff
    
    mcf_filtered <- mcf_filtered %>%
      rownames_to_column(var="Bin") %>%
      as_tibble() %>% #NB - sometimes the select across stuff only works on tibbles and not df...
      select(where(~ any(.x >= cutoff))) %>%  #select(where(~!all(sum(.x)==0))
      column_to_rownames(var="Bin")
    
    
    if (length(colnames(mcf_filtered)) == 0) {
      ggplot() + annotate("text", x=4, y=25, size=8, label="Not enough data to generate a plot with this module filter.") + theme_void()
    
    } else {
  
      #Make pretty text for modules - taken from the METADATA
      KO_short_names <- as.data.frame(cbind(mcf_data$raw_results$METADATA[,1], mcf_data$raw_results$METADATA[,3]))
      colnames(KO_short_names) <- c("ID", "Short name")
      KO_short_names$Name <- apply(KO_short_names[, 1:2], 1, paste, collapse = " - ")
      
      #Add these pretty names back to the matrix columns
      mcf_names <- as.data.frame(colnames(mcf_filtered)) %>%
        left_join(KO_short_names, by = c("colnames(mcf_filtered)" = "ID")) %>%
        select(Name)
      colnames(mcf_filtered) <- mcf_names$Name
      
      
      ord_x <- 1:length(colnames(mcf_filtered))
      if (length(colnames(mcf_filtered)) > 2) {
        ord_x <- hclust(dist(scale(t(mcf_filtered)), method = "euclidean"), method = "ward.D")$order
      }
      
  
          #Generate heatmap - transposed version
      # plot <- pheatmap(as.data.frame(mcf_filtered),
      #                  cellwidth = 10,
      #                  cellheight = 8,
      #                  fontsize = 6,
      #                  cluster_rows = cluster_rows,
      #                  cluster_cols = cluster_cols,
      #                  filename = "mcf.pdf")
      
      write.table(mcf_filtered, file="mcf.tsv", sep="\t", quote=FALSE)
      
          
      mcf_filtered_ext <- mcf_filtered %>% rownames_to_column(var="Bin") %>% pivot_longer(!Bin, names_to = "Name", values_to = "Value")
      
      
      plot <- ggplot(mcf_filtered_ext) +
        geom_tile(aes(y=Bin, x=Name, fill=Value), color="white") +
        labs(title = "Module completion fraction", x="Module name", y="MAGs") +
        theme(legend.position = 'top',
              axis.text.x = element_text(size = 10, angle = 90, vjust = 1, hjust=1),
              axis.text.y = element_text(size = 10)) +
        guides(fill = guide_colourbar(title="", barwidth = 20, barheight = 0.5)) +
        scale_y_discrete(position = "right") +
        scale_x_discrete(limits = mcf_filtered_ext$Name[ord_x])
  
  
      ggsave("mcf.pdf", plot)
      plot
    }
  })
  
  output$downloadMCFPlot <- downloadHandler(
    filename = function() {"mcf.pdf"},
    content = function(file) {file.copy("mcf.pdf", file, overwrite=TRUE)}
  )
  
  output$downloadMCFPlot_data <- downloadHandler(
    filename = function() {"mcf.tsv"},
    content = function(file) {file.copy("mcf.tsv", file, overwrite=TRUE)}
  )
  
  

}


shinyApp(ui = ui_main, server = server)
