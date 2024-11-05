# List of required packages
required_packages <- c(
  "shiny", 
  "shinydashboard", 
  "shinydashboardPlus", 
  "shinycssloaders", 
  "biomaRt", 
  "rtracklayer", 
  "BSgenome.Hsapiens.UCSC.hg38", 
  "Gviz"
)

# Function to check and install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Install Bioconductor packages separately if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("rtracklayer", update = FALSE, ask = FALSE)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update = FALSE, ask = FALSE)
BiocManager::install("Gviz", update = FALSE, ask = FALSE)

# Install CRAN packages
install_if_missing(required_packages)

# Load all packages
lapply(required_packages, library, character.only = TRUE)

cat("Setup complete. All required packages are installed and loaded.\n")
