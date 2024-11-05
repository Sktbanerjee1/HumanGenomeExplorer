library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinycssloaders)

# Define the UI for the app using shinydashboard and shinydashboardPlus
ui <- dashboardPage(
  title = "HGE",
  header = dashboardHeader(
    title = "HGE",
    titleWidth = 250
  ),
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("Human Genome", tabName = "genome", icon = icon("dna")),
      menuItem("Region View", tabName = "detailed", icon = icon("search"))
    )
  ),
  body = dashboardBody(
    tabItems(
      # Tab for displaying all chromosome ideograms
      tabItem(
        tabName = "genome",
        fluidRow(
          box(
            title = "History of Human Genome Sequencing",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            p("The history of human genome sequencing spans several decades of scientific advancement:"),
            tags$ul(
              tags$li(
                strong("1984:"), 
                " The idea of sequencing the entire human genome is proposed during the Santa Fe workshop."
              ),
              tags$li(
                strong("1990:"), 
                " The Human Genome Project (HGP) officially begins, with the goal of sequencing the entire human genome within 15 years."
              ),
              tags$li(
                strong("2001:"), 
                " A draft sequence of the human genome is published by the HGP and Celera Genomics, marking a major milestone."
              ),
              tags$li(
                strong("2003:"), 
                " The Human Genome Project is completed, providing a nearly complete sequence of the human genome, covering approximately 99% of the genome with 99.99% accuracy."
              ),
              tags$li(
                strong("2010:"), 
                " The 1000 Genomes Project is launched, aiming to sequence the genomes of a large number of individuals to study genetic variation."
              ),
              tags$li(
                strong("2022:"), 
                " The Telomere-to-Telomere (T2T) Consortium publishes the first truly complete sequence of a human genome, including previously missing regions."
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "All Human Chromosomes",
            width = 12,
            status = "primary",
            withSpinner(plotOutput("allChromosomesPlot", height = "800px"))
          ),
          box(
            title = "Ideogram Track Description",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            p("The ideogram track provides a graphical representation of each chromosome."),
            p("Dark and light bands indicate regions of different staining intensities with Giemsa stain:"),
            tags$ul(
              tags$li("Dark Bands: Generally correspond to regions of lower gene density and are rich in A-T base pairs."),
              tags$li("Light Bands: Tend to be gene-rich regions and are GC-rich.")
            ),
            p("The centromere, represented as a constriction in the ideogram, is the region where the two chromatids are joined.")
          )
        )
      ),
      # Tab for detailed view of selected chromosome
      tabItem(
        tabName = "detailed",
        fluidRow(
          box(
            title = "Selection Panel",
            width = 4,
            status = "info",
            selectInput(
              inputId = "popularGene",
              label = "Jump to Popular Gene:",
              choices = c("None", "TP53", "BRCA1", "EGFR", "MYC"),
              selected = "None"
            ),
            textInput("geneSearch", "Search Gene Name:"),
            actionButton("searchGene", "Search"),
            selectInput(
              inputId = "chromosome",
              label = "Select Chromosome:",
              choices = paste0("chr", c(1:22, "X", "Y")),
              selected = "chr1"
            ),
            numericInput(
              inputId = "start",
              label = "Selection Start (in base pairs):",
              value = 1,
              min = 1
            ),
            numericInput(
              inputId = "end",
              label = "Selection End (in base pairs):",
              value = 1000000,
              min = 1
            ),
            actionButton("zoomIn", "Zoom In"),
            actionButton("zoomOut", "Zoom Out"),
            actionButton("update", "Update Plot"),
            textOutput("warningMessage")
          ),
          box(
            title = "Detailed Genome Plot",
            width = 8,
            status = "primary",
            withSpinner(plotOutput("detailedPlot", height = "800px"))
          )
        ),
        fluidRow(
          box(
            title = "Genome Features Description",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            p("This detailed view provides information on three main features of the genome:"),
            tags$ul(
              tags$li(
                strong("Gene Track:"), 
                " Displays annotated genes and their structures, including exons and introns."
              ),
              tags$li(
                strong("SNP Track:"), 
                " Shows Single Nucleotide Polymorphisms (SNPs) within the selected region, with allele identifiers."
              ),
              tags$li(
                strong("GC Content Track:"), 
                " Illustrates the GC content in the selected genomic region, providing insight into the base composition."
              )
            )
          )
        )
      )
    )
  ),
  footer = dashboardFooter(left = "Human Genome Explorer © 2024")
)
