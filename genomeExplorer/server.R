library(shiny)
library(Gviz)
library(biomaRt)
library(gridExtra)
library(grid)
library(BSgenome.Hsapiens.UCSC.hg38)
library(trackViewer)
library(rtracklayer)

# Load the modules
source("modules/allChromosomesModule.R")
source("modules/detailedViewModule.R")

# Define the main server logic
server <- function(input, output, session) {
  # Set up biomaRt to fetch gene information
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  # Set up a new BioMart object to fetch SNP information
  human_variation <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
  # Reactive values for shared data across modules
  plot_params <- reactiveValues(chromosome = "chr1", start = 1, end = 1e6)
  
  # Call the modules
  allChromosomesModule(input, output, session)
  detailedViewModule(input, output, session, plot_params, ensembl, human_variation)
}
