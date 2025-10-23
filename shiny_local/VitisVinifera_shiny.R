library(shiny)
library(shinythemes)
library(dplyr)
library(DT)
library(ggplot2)
library(shinyjs)
library(shinyWidgets)
library(data.table)
library(shinycssloaders)
library(writexl)
library(Biostrings)

metadata <- readxl::read_xlsx("26169genes_with_AthalHomologs_allIDs_exprPatterns.xlsx")
Rlogs <- read.table("Rlogs.csv", sep="\t", header = TRUE)
rawCounts <- read.table("RawCounts.csv", sep="\t", header = TRUE)

ui <- fluidPage(
  useShinyjs(),
    tags$div(
      class = "shiny-title-panel",
      tags$h1(
        "Grapevine Guardians: Unraveling Rpv Pyramidization’s Impact on Immunity",
        style = paste(
          "color: #2E8B57;",
          "font-family: 'Segoe UI', Tahoma, sans-serif;",
          "font-weight: 300;",
          "font-size: 2.5rem;"
        )
      )
    ),
  tags$br(),
  tags$br(),
  fluidRow(
    column(width = 1, align = "center", icon("wine-glass", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-bottle", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-glass", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-bottle", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-glass", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-bottle", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-glass", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-bottle", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-glass", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-bottle", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-glass", class = "fa-3x")),
    column(width = 1, align = "center", icon("wine-bottle", class = "fa-3x"))
  ),
  tags$br(),
  tags$br(),
  mainPanel(
    tabsetPanel(id = "plateFill",
                tabPanel("Data overview", icon=icon(name="magnifying-glass", lib="font-awesome"), value="otherTab",
                         DTOutput("conversionTab"), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(),tags$br()),
                
                tabPanel( title = "Microview", icon=icon(name="magnifying-glass-plus", lib="font-awesome"), value="detailsTab",
                          tags$br(),
                          tableOutput("rowDetails"),
                          tags$br(),
                          fluidRow(
                            column(6, align="right",  downloadButton("downloadCDSseq", "Download CDS (nt, fasta)")),
                            column(6, align="left",   downloadButton("downloadProtseq", "Download Protein (aa, fasta)"))
                          ),
                          tags$br(),
                          tags$br(),
                          fluidRow(
                            column(4, align="right",  uiOutput("UniprotLink") ),
                            column(4, align="centre", uiOutput("Ncbiacc") ),
                            column(4, align="left",   uiOutput("TAIR10protLink") )
                          ),
                          tags$br(),
                          tags$br(),
                          tags$h3("Data analysis"),
                          tags$hr(style = "border-top: 1px solid dimgray;"),
                          tags$br(),
                          tableOutput("rlogsDetail"),
                          tags$br(),
                          plotOutput("BoxPlotRlogs"),
                          tags$br(),
                          tags$br(),
                          fluidRow(
                            column(width=12, align="left", div(
                              div(style="display: inline-block; width: 425px ;", tags$img(src="legend.png",
                                                                                          width = "400px",
                                                                                          style="vertical-align: top;")),
                              
                            )),
                            tags$br(),
                            tags$br(),
                            "The PCA plot of a co-transcriptional module is shown only if the focal protein-coding gene is significantly differentially expressed (|log₂FC| > 1, FDR < 0.05) in at least one condition and if it correlates (Person’s ρ > 0.817) with more than two other protein-coding genes.",
                            tags$br(),
                            tags$br(),
                            column(width=12, align="centre", 
                            div(style="display: inline-block; width: 425px ;", imageOutput("VitisID_plot", width = "400px") ))),
                          
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          fluidRow(
                            column(width=12, align="left",
                          "If a co-transcriptional module is identified, we calculate the Spearman correlation coefficient between the first principal component (derived from rlog-transformed transcriptional values across four genotypes, three biological replicates, and three time points) and the simplified expression pattern of the focal gene (coded as 0 for no significant change compared to the reference, −1 for significant downregulation, and +1 for significant upregulation). P-values of all computed Spearman correlation coefficients were retained and adjusted for multiple testing (FDR, Bonferroni). The adjusted values are reported in the table below.",
                          tags$br(),
                          tags$br(),
                          tableOutput("moduleInfo_VitisID"),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          "List of cotranscriptional partners:",
                          tags$br(),
                          tags$br(),
                          DTOutput("partners_VitisID"),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          column(width=12, align="centre", downloadButton(outputId = "downloadCotrP", label="Cotranscriptional partners (csv)")),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br(),
                          tags$br()
                          )),
                          
                      
                ),
                
                tabPanel("Abbreviations", icon=icon(name="star-of-life", lib="font-awesome"),
                         tags$br(),
                         tags$span("Used color code:", style="font-size:18px; color:dimgray; font-weight:bold"),
                         tags$br(),
                         div(style="display: inline-block; width: 800px ;", tags$img(src="legend.png",
                                                                                     width="750px",
                                                                                     style="vertical-align: middle;") ),
                         tags$br(),
                        
                         tags$span("the Rpv12 genotype - goldenrod", style="font-size:18px; color:goldenrod; font-weight:normal"),
                         tags$br(),
                         tags$span("the Rpv12+1 genotype - salmon", style="font-size:18px; color:salmon; font-weight:normal"),
                         tags$br(),
                         tags$span("the Rpv12+1+3 genotype - cornflowerblue", style="font-size:18px; color:cornflowerblue; font-weight:normal"),
                         tags$br(),
                         tags$span("the susceptible genotype (Pinot noir, the control, the reference) - dimgray", style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br(),
                         tags$br(),
                         tags$span("Symbols for the time points in plots:", style="font-size:18px; color:dimgray; font-weight:bold"),
                         tags$br(),
                         tags$div(
                           tags$span(class = "circle"), " R pch symbol 20 - Circle represents data at 0 hours post infection (hpi).", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         # Triangle symbol and description
                         tags$div(
                           tags$span(class = "triangle"), " R pch symbol 17 - Triangle represents data at 6 hpi.", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         # Square symbol and description
                         tags$div(
                          tags$span(class = "square"), " R pch symbol 15 - Square represents data at 24 hpi.", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span("Symbols for the transcription pattern:", style="font-size:18px; color:dimgray; font-weight:bold"),
                         tags$br(),
                         tags$span(
                           "E - EUQL - equal to the control genotype at the given time", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "U - UP - upregulated - cut-off: FDR<0.05, log2FCH>1", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "D - DOWN - downregulated - cut-off: FDR<0.05, log2FCH<-1", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$br(),
                         tags$span("Transcription patterns are encoded into 9 letters long codes:", style="font-size:18px; color:dimgray; font-weight:bold"),
                         tags$br(),
                         tags$span(
                           "Rpv12 (0 hpi, 6hpi, 24 hpi) | Rpv12+1 (0 hpi, 6hpi, 24 hpi) | Rpv12+1+3 (0 hpi, 6hpi, 24 hpi)", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "for example: DDD EEE DDE", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$br(),
                         tags$span("Gene categorization:", style="font-size:18px; color:dimgray; font-weight:bold"),
                         tags$br(),
                         tags$span(
                           "IEV - Initial Expression Variations - changes detected at 0 hpi or 0 and 6 hpi (transcription pattern: UUE, UEE, DDE, DEE)", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "ER - Early Response - changes detected at 6 hpi and 24 hpi (transcription pattern: EUU, EDD)", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "TRS - Transient Response on Stress - changes detected at 6 hpi and 24 hpi (transcription pattern: EUU, EDD)", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "LT - Late Response - changes detected at 24 hpi (transcription pattern: EEU, EED)", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "SCh - Sustained Change - changes detected in all three time points (transcription pattern: UUU, DDD)", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "CPs - Complex Patterns - changes different from the abovementioned (examples: UEU, DED, DEU, DDU)", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$br(),
                         tags$span("Grouping of gene categories:", style="font-size:18px; color:dimgray; font-weight:bold"),
                         tags$br(),
                         tags$span(
                           "Group I - Rpv12-Triggered Group - a gene category shared by all three cultivars", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "Group II - Rpv3-Attenuated Group - a gene category shared by Rpv12 and Rpv12+1 cultivars", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "Group III - Rpv1-Triggered Group - a gene category shared by all three cultivars", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "Group IV - a recovered gene category - patterns shared by Rpv12 and Rpv12+1+3 cultivars", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span(
                           "Group V - cultivar-specific gene categories (changes found just in one cultivar)", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$span("Group V has three subgroups: (a) Rpv12-specific changes; (b) Rpv12+1-specific changes; (c) Rpv12+1+3-specific changes", style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br(),
                         tags$span(
                           "Group VI - complex patterns across cultivars ", style="font-size:18px; color:dimgray; font-weight:normal"
                         ),
                         tags$br(),
                         tags$img(src = "gene_category.png", width = "800px", height = "200px"),
                         tags$br(),
                         tags$br(),
                         tags$span("Sample IDs:", style="font-size:18px; color:dimgray; font-weight:bold"),
                         tags$br(),
                         tags$span("A, B, C - replicates", style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br(),
                         tags$span(".0 - 0 hpi", style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br(),
                         tags$span(".6 - 6 hpi", style="font-size:18px; color:dimgray; font-weight:normal"),  
                         tags$br(),
                         tags$span(".24 - 24 hpi", style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br(),
                         tags$span("example: sample ID Rpv12+1.0.A means the genotype RVP12+1 collected at 0 hpi, replicate A", style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br(),
                         tags$br(),
                         tags$span("UniProt IDs", style="font-size:18px; color:dimgray; font-weight:bold"),
                         tags$br(),
                         tags$span("are compatible with following databases to perform gene ontology analysis (GOA), statistical overrepresentation tests, or just functional classification: ", 
                                   style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br(),
                         tags$a(href ="https://string-db.org/", ">> STRING-DB", 
                                style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br(),
                         tags$span(" and ", style="font-size:10px; color:dimgray"),
                         tags$br(),
                         tags$a(href ="https://www.pantherdb.org/", ">> PANTHER-DB", 
                                style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br(),
                         tags$br(),
                         tags$br(),
                         tags$br(),
                         tags$br(),
                         tags$br(),
                         tags$br()
                           ),
                
                tabPanel("Bulk downloads", icon=icon(name="download", lib="font-awesome"),
                 tags$br(),
                 fluidRow(
                  tags$span("Data overview table (information about 26 169 detected transcripts):", style="font-size:18px; color:dimgray; font-weight:bold"),
                  tags$br(),
                  column(downloadButton(outputId = "downloadMetadata", label = "Metadata"), width=12),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$span("Matrix with normalized values (DESeq2 rlog function and SVA parametric combat function):", style="font-size:18px; color:dimgray; font-weight:bold"),
                  tags$br(),
                  column(downloadButton(outputId = "downloadRlogs", label="rlog values"), width=12),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$span("Matrix with raw counts (Subread featureCount function):", style="font-size:18px; color:dimgray; font-weight:bold"),
                  tags$br(),
                  column(downloadButton(outputId = "downloadRawCounts", label="Raw counts"), width=12),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$span("Reference fasta file can be found here:", style="font-size:18px; color:dimgray; font-weight:bold"),
                  tags$a(href ="https://grapedia.org/files-download/", "PN40024.v4", 
                         style="font-size:18px; color:dimgray; font-weight:normal"),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$br(),
                  tags$br()
                 )        
                           ),
                tabPanel("Data analysis", icon=icon(name="github", lib="font-awesome"),
                         fluidRow(
                         tags$a(href ="https://github.com/vierocka/plant_immunity_Vitis", ">> GitHub link", 
                                style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br()
                         )        
                ),
                tabPanel("Project information", icon = icon(name="circle-info", lib="font-awesome"),
                         fluidRow(
                         tags$h4("Acknowledgement"),
                         tags$br(),
                         tags$br(),
                        
                         fluidRow(column(width=12, align="center", div(
                           div(style="display: inline-block; width: 150px ;", tags$img(src="Vinselekt_Michlovsky.png",
                                                                                       height="75px",
                                                                                       style="vertical-align: middle;") ),
                           
                           div(style="display: inline-block; width: 150px ;", tags$img(src="Mendelu_logo_2019.jpg",
                                                                                      height="75px",
                                                                                      style="vertical-align: middle;")),
                           
                           div(style="display: inline-block; width: 150px ;", tags$img(src="UoC.png",
                                                                                      height="75px",
                                                                                      style="vertical-align: middle;"))
                           
                         ))),
                         tags$br(),
                         tags$br(),
                         span("Contact:", style="font-size:18px; color:dimgray; font-weight:bold"),
                         tags$br(),
                         tags$a(href ="https://www.linkedin.com/in/viera-kovacova-35927877/", ">> Viera Kovacova, PhD", 
                                style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br(),
                         tags$a(href ="https://is.mendelu.cz/lide/clovek.pl?id=3560", ">> Prof. Miroslav Baranek, PhD", 
                                style="font-size:18px; color:dimgray; font-weight:normal"),
                         tags$br()
                                )
                         )
                         )
                ),

# Back to top button
tags$style(HTML("
    #backToTop {
      position: fixed;
      bottom: 20px;
      right: 20px;
      z-index: 1000;
      display: none;  /* hidden by default */
    }
  ")),

actionButton("backToTop", "⬆ Back to Top"),

# JavaScript to show button when scrolling down and scroll to top on click
tags$script(HTML("
    $(document).ready(function(){
      // Show or hide button
      $(window).scroll(function(){
        if ($(this).scrollTop() > 100) {
          $('#backToTop').fadeIn();
        } else {
          $('#backToTop').fadeOut();
        }
      });
      
      // Scroll to top on click
      $('#backToTop').click(function(){
        $('html, body').animate({scrollTop : 0}, 400);
        return false;
      });
    });
  "))
    )
  

server <- function(input, output, session) {
  Sys.sleep(1) 
  
  output$conversionTab <- renderDT({ metadata })
  
  output$downloadMetadata <- downloadHandler(
    filename = function() { "metadata.csv" },
    content = function(file) {
      write.csv(metadata, file, row.names = FALSE)
    }
  )
  

  output$downloadRlogs <- downloadHandler(
    filename = function() { "data_files/Rlogs.csv" },
    content = function(file) {
      write.csv(Rlogs, file, row.names = FALSE)
    }
  )
  
  output$downloadRawCounts <- downloadHandler(
    filename = function() { "data_files/rawCounts.csv" },
    content = function(file) {
      write.csv(rawCounts, file, row.names = FALSE)
    }
  )
  

  # Hide the "Conditional" tab on startup
  hideTab(inputId = "plateFill", target = "detailsTab")
  
# condition for zooming into data after row click
  observeEvent(input$conversionTab_rows_selected, {
    req(input$conversionTab_rows_selected)
    selected_row <- input$conversionTab_rows_selected
    if (selected_row != ""){
    selected_data <- metadata[selected_row,]
    VitisID <- metadata[selected_row,1]
    UniProt <- metadata[selected_row,2]
    NCBIacc <-  metadata[selected_row,8]
    Tair10 <-  metadata[selected_row,16]
    
    output$rowDetails <- renderTable({ selected_data })

    output$UniprotLink <- renderUI({
      tags$a(href = paste0("https://www.uniprot.org/uniprotkb/", UniProt, "/entry"),
             "UniProt link",
             style="font-size:18px; color:dimgray; font-weight:normal",
             target = "_blank")
    })
    
    output$Ncbiacc <- renderUI({
      tags$a(href = paste0("https://www.ncbi.nlm.nih.gov/protein/", NCBIacc),
             "NCBI ASM3070453v1 protein accession",
             style="font-size:18px; color:dimgray; font-weight:normal",
             target = "_blank")
    })
    
    output$TAIR10protLink <- renderUI({
      tags$a(href = paste0("https://www.arabidopsis.org/results?mainType=general&searchText=", Tair10, "&category=genes"),
             "TAIR10 link",
             style="font-size:18px; color:dimgray; font-weight:normal",
             target = "_blank")
    })
    
    output$rlogsDetail <- renderTable({ Rlogs[match(VitisID, Rlogs$geneID),] })
    
    output$BoxPlotRlogs <- renderPlot({
    BoxPlData <- matrix(unlist(c(Rlogs[match(VitisID, Rlogs$geneID), as.double(c(2,14,26))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(3,15,27))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(4,16,28))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(5,17,29))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(6,18,30))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(7,19,31))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(8,20,32))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(9,21,33))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(10,22,34))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(11,23,35))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(12,24,36))], Rlogs[match(VitisID, Rlogs$geneID), as.double(c(13,25,37))])), ncol=12)
    colnames(BoxPlData) <- c("Rpv12_0hpi","Rpv12_6hpi","Rpv12_24hpi","Rpv12+1_0hpi","Rpv12+1_6hpi","Rpv12+1_24hpi", "Rpv12+1+3_0hpi","Rpv12+1+3_6hpi","Rpv12+1+3_24hpi", "susceptible_0hpi", "susceptible_6hpi", "susceptible_24hpi")
    par(mar=c(8,3.5,1.5,1), mgp=c(2.5,0.75, 0), cex.lab=1, cex.axis=1, cex.main=0.9)
    boxplot(BoxPlData, main=VitisID, xaxt="n", xlab="", ylab="rlog values", col=c("goldenrod","goldenrod","goldenrod","salmon","salmon","salmon","cornflowerblue","cornflowerblue","cornflowerblue","dimgray","dimgray","dimgray"))
    axis(1, at=c(1:12), labels = c("Rpv12_0hpi","Rpv12_6hpi","Rpv12_24hpi","Rpv12+1_0hpi","Rpv12+1_6hpi","Rpv12+1_24hpi", "Rpv12+1+3_0hpi","Rpv12+1+3_6hpi","Rpv12+1+3_24hpi", "susceptible_0hpi", "susceptible_6hpi", "susceptible_24hpi"), las=2)
    points(x=c(1:12), y=BoxPlData[1,], pch=rep(c(20,17,15), 4), cex=1.25, col="gray60")
    points(x=c(1:12), y=BoxPlData[2,], pch=rep(c(20,17,15), 4), cex=1.25, col="gray60")
    points(x=c(1:12), y=BoxPlData[3,], pch=rep(c(20,17,15), 4), cex=1.25, col="gray60") })
    
    CDSes <- readDNAStringSet("Vitis_vinifera.PN40024.v4.cds.all.fa")
    Proteins <- readAAStringSet("Vitis_vinifera.PN40024.v4.pep.all.fa")
    
    output$downloadCDSseq <- downloadHandler(
      filename = function() { paste0(VitisID, "_t001_CDS_Nt.fa") },
      content = function(file) {
        writeXStringSet( CDSes[ grep( paste0(VitisID, "_t001"), names(CDSes)) ], file )
      }
    )
    
    output$downloadProtseq <- downloadHandler(
      filename = function() { paste0(VitisID, "_P001_AA.fa") },
      content = function(file) {
        writeXStringSet( Proteins[ grep( paste0(VitisID, "_P001"), names(Proteins)) ], file )
      }
    )

    if (file.exists( paste0("www/", VitisID, "_module_PCA_CorrCutOff08.svg") )) {
    output$VitisID_plot <- renderImage({
       list(
        src = paste0("www/", VitisID, "_module_PCA_CorrCutOff08.svg"),
        contentType = "image/svg+xml",
        width = 500
      )
    }, deleteFile = FALSE) }else{
      
      output$VitisID_plot <- renderImage({
        list(
          src = paste0("www/", "NoPCA.svg"),
          contentType = "image/svg+xml",
          width = 500
        )
      }, deleteFile = FALSE)
      
    }
    
    moduleInfo <- read.table("TranscrpModule_Spearman_Pval_Padj_FDR.csv", sep=",", header = TRUE)
    output$moduleInfo_VitisID <-  renderTable({ moduleInfo[match(VitisID, moduleInfo$geneID),] }, digits = 6)
    if (file.exists(paste("CotranscModules/", VitisID, "_info.csv", sep = ""))) {
    CotrPart <- read.table(paste("CotranscModules/", VitisID, "_info.csv", sep = ""), sep="\t",  header = TRUE)
    MatchAbbrev <- as.data.frame( cbind(CotrPart, metadata[match(rownames(CotrPart), metadata$PN40024v4_ENSMBL),c(2,3,7,8,9)]))}else{
      MatchAbbrev <- c()
    }
    output$partners_VitisID <-  renderDT({ MatchAbbrev })
    
    output$downloadCotrP <- downloadHandler(
      filename = function() { paste(VitisID, "_cotranscrPartners.csv", sep = "") },
      content = function(file) {
        write.csv(MatchAbbrev, file, row.names = TRUE )
      }
    )
    
    showTab(inputId = "plateFill", target = "detailsTab")
    updateTabsetPanel(session, "plateFill", selected = "detailsTab")
    
    }else{
      updateTabsetPanel(session, "plateFill", selected = "otherTab") 
      hideTab(inputId = "plateFill", target = "detailsTab")
    }
  })
  
  table_proxy <- dataTableProxy("conversionTab")
  # When the user navigates tabs
  observeEvent(input$plateFill, {
    if ( !identical(input$plateFill, "detailsTab")){
      # Hide the Conditional tab again
      hideTab(inputId = "plateFill", target = "detailsTab")
      
      # Clear the detail content and reset selection
      output$selected_data <- renderTable(NULL)
      selectRows(table_proxy, NULL)
    }
  })

}

setwd("~/Dropbox/MendelUni_Vinselect/shiny_local/")
# library(rsconnect)
# deployApp(appName = "playing_with_immunity", appDir = "~/Dropbox/MendelUni_Vinselect/shinyData/")
shinyApp(ui = ui, server = server)


