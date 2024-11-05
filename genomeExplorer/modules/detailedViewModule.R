detailedViewModule <- function(input, output, session, plot_params, ensembl, human_variation) {
  # Increase the timeout limit for biomaRt queries
  options(timeout = 120)
  
  # Pre-computed coordinates for popular genes
  precomputed_genes <- list(
    "TP53" = list(chromosome = "chr17", start = 7668402 - 10000, end = 7668402 + 10000),
    "BRCA1" = list(chromosome = "chr17", start = 43044295 - 10000, end = 43044295 + 10000),
    "EGFR" = list(chromosome = "chr7", start = 55086714 - 10000, end = 55086714 + 10000),
    "MYC" = list(chromosome = "chr8", start = 127735434 - 10000, end = 127735434 + 10000)
  )
  
  # Update for popular genes using pre-computed coordinates
  observeEvent(input$popularGene, {
    if (input$popularGene != "None") {
      gene_info <- precomputed_genes[[input$popularGene]]
      plot_params$chromosome <- gene_info$chromosome
      plot_params$start <- gene_info$start
      plot_params$end <- gene_info$end
      updateSelectInput(session, "chromosome", selected = gene_info$chromosome)
      updateNumericInput(session, "start", value = gene_info$start)
      updateNumericInput(session, "end", value = gene_info$end)
    }
  })
  
  # Gene search functionality with a timeout handler
  observeEvent(input$searchGene, {
    gene_name <- input$geneSearch
    if (gene_name != "") {
      tryCatch({
        gene_data <- getBM(
          attributes = c("chromosome_name", "start_position", "end_position"),
          filters = "hgnc_symbol",
          values = gene_name,
          mart = ensembl
        )
        
        if (nrow(gene_data) > 0) {
          chromosome <- paste0("chr", gene_data$chromosome_name[1])
          start_pos <- gene_data$start_position[1]
          end_pos <- gene_data$end_position[1]
          center <- (start_pos + end_pos) / 2
          zoom_size <- 20000  # Set the zoom size to 20,000 bp
          new_start <- max(1, round(center - zoom_size / 2))
          new_end <- round(center + zoom_size / 2)
          
          plot_params$chromosome <- chromosome
          plot_params$start <- new_start
          plot_params$end <- new_end
          updateSelectInput(session, "chromosome", selected = chromosome)
          updateNumericInput(session, "start", value = new_start)
          updateNumericInput(session, "end", value = new_end)
        } else {
          showNotification("Gene not found. Please check the gene name and try again.", type = "error")
        }
      }, error = function(e) {
        showNotification("Error fetching gene data. Please try again later.", type = "error")
      })
    }
  })
  
  # Update plot functionality with minimum region size constraint
  observeEvent(input$update, {
    plot_params$chromosome <- input$chromosome
    plot_params$start <- max(1, input$start)  # Ensure start is at least 1
    plot_params$end <- max(input$start + 1000, input$end)  # Ensure end is at least 1,000 bp greater than start
    updateNumericInput(session, "start", value = plot_params$start)
    updateNumericInput(session, "end", value = plot_params$end)
  })
  
  # Ensure zoom in and zoom out functionalities respect the minimum region size
  observeEvent(input$zoomIn, {
    region_size <- plot_params$end - plot_params$start
    zoom_factor <- 0.5  # Zoom in by reducing the region size by half
    new_size <- max(1000, round(region_size * zoom_factor))  # Minimum size of 1,000 bp
    center <- (plot_params$start + plot_params$end) / 2
    new_start <- max(1, round(center - new_size / 2))
    new_end <- round(center + new_size / 2)
    
    plot_params$start <- new_start
    plot_params$end <- new_end
    updateNumericInput(session, "start", value = new_start)
    updateNumericInput(session, "end", value = new_end)
  })
  
  observeEvent(input$zoomOut, {
    region_size <- plot_params$end - plot_params$start
    zoom_factor <- 2  # Zoom out by doubling the region size
    new_size <- min(5e6, round(region_size * zoom_factor))  # Maximum size of 5,000,000 bp
    new_size <- max(1000, new_size)  # Ensure the new size is at least 1,000 bp
    center <- (plot_params$start + plot_params$end) / 2
    new_start <- max(1, round(center - new_size / 2))
    new_end <- round(center + new_size / 2)
    
    plot_params$start <- new_start
    plot_params$end <- new_end
    updateNumericInput(session, "start", value = new_start)
    updateNumericInput(session, "end", value = new_end)
  })
  
  # Compute GC content
  compute_gc_content <- reactive({
    if (is.null(plot_params$chromosome)) return(NULL)
    
    # Get the sequence for the specified region
    seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                  names = plot_params$chromosome, 
                  start = plot_params$start, 
                  end = plot_params$end)
    
    # Calculate GC content using a sliding window
    window_size <- 100  # Size of the sliding window
    gc_content <- letterFrequencyInSlidingView(seq, view.width = window_size, letters = "GC", as.prob = TRUE)
    
    # Adjust the lengths to match
    positions <- seq(plot_params$start + window_size / 2, plot_params$end - window_size / 2, by = 1)
    if (length(positions) > length(gc_content)) {
      positions <- positions[seq_along(gc_content)]
    } else if (length(gc_content) > length(positions)) {
      gc_content <- gc_content[seq_along(positions)]
    }
    
    # Create a data frame for GC content
    gc_df <- data.frame(position = positions, gc_content = gc_content)
    return(gc_df)
  })
  
  # Fetch SNP data using biomaRt with a timeout handler
  fetch_snp_data <- reactive({
    if (is.null(plot_params$chromosome)) return(NULL)
    
    region_size <- plot_params$end - plot_params$start
    if (region_size > 1e5) return(NULL)  # Skip fetching SNPs if the region is larger than 1 Mb
    
    tryCatch({
      snp_data <- getBM(
        attributes = c("refsnp_id", "chr_name", "chrom_start"),
        filters = c("chr_name", "start", "end"),
        values = list(sub("chr", "", plot_params$chromosome), plot_params$start, plot_params$end),
        mart = human_variation
      )
      
      if (nrow(snp_data) > 0) {
        snp_data$chrom_start <- as.numeric(snp_data$chrom_start)
        return(snp_data)
      }
      return(NULL)
    }, error = function(e) {
      showNotification("Error fetching SNP data. Please try again later.", type = "error")
      return(NULL)
    })
  })
  
  # Render the detailed plot
  output$detailedPlot <- renderPlot({
    # Error handling for large region size
    validate(
      need(plot_params$start < plot_params$end, "Invalid selection range. Please adjust the start and end positions."),
      need(plot_params$end - plot_params$start <= 5e6, "Selected region is too large. Please select a smaller region (less than 5 Mb).")
    )
    
    region_size <- plot_params$end - plot_params$start
    
    # Create a GenomeAxisTrack to display positions in megabases
    genome_axis_track <- GenomeAxisTrack(
      labelPos = "above",
      fontsize = 18,
      name = "Position",
      height = 0.5  # Adjusted height
    )
    
    # Create an IdeogramTrack for the selected chromosome
    ideogram_track <- IdeogramTrack(
      genome = "hg38",
      chromosome = plot_params$chromosome,
      showBandId = FALSE,
      showId = FALSE,
      height = 0.7  # Adjusted height
    )
    
    # Create a BiomartGeneRegionTrack with ENSEMBL annotations
    gene_track <- BiomartGeneRegionTrack(
      genome = "hg38",
      chromosome = sub("chr", "", plot_params$chromosome),
      start = plot_params$start,
      end = plot_params$end,
      name = "ENSEMBL",
      biomart = ensembl,
      transcriptAnnotation = "symbol",
      collapseTranscripts = "gene",
      fill = "orange",
      col = "black",
      col.line = "darkgray",
      lwd = 1,
      cex = 1.5,
      fontcolor.group = "black",
      background.title = "black",
      height = 1  # Adjusted height
    )
    
    # Prepare plot tracks
    plot_tracks <- list(genome_axis_track, ideogram_track, gene_track)
    
    # Add SNP track if the region is not too large
    if (region_size <= 1e6) {
      snp_data <- fetch_snp_data()
      
      if (!is.null(snp_data) && nrow(snp_data) > 0) {
        # Optionally, limit the number of SNPs for visualization
        if (nrow(snp_data) > 1000) {  # Example threshold: 1,000 SNPs
          snp_data <- snp_data[seq(1, nrow(snp_data), length.out = 1000), ]
        }
        
        snp_track <- AnnotationTrack(
          start = snp_data$chrom_start,
          end = snp_data$chrom_start,
          chromosome = plot_params$chromosome,
          genome = "hg38",
          id = snp_data$refsnp_id,
          name = "SNPs",
          col = "black",
          fill = "darkred",
          background.title = "black",
          stacking = "dense",  # Use "dense" to pack features tightly
          redundant = TRUE,  # Allow redundancy to avoid excess stacks
          group = snp_data$refsnp_id,  # Group SNPs by their IDs
          groupAnnotation = "id",  # Display allele IDs for SNPs
          height = 0.3  # Reduced height
        )
        
        plot_tracks <- c(plot_tracks, list(snp_track))
      }
    } else {
      showNotification("SNP track can't be plotted because the region is too large. Please zoom in.", type = "warning")
    }
    
    # Add GC content track if data is available
    gc_df <- compute_gc_content()
    if (!is.null(gc_df) && nrow(gc_df) > 0) {
      gc_track <- DataTrack(
        range = IRanges(start = gc_df$position, width = 1),
        data = gc_df$gc_content,
        genome = "hg38",
        chromosome = plot_params$chromosome,
        name = "GC Content",
        type = "h",  # Histogram
        col.histogram = "blue",
        fill.histogram = "lightblue",
        background.title = "black",
        height = 1  # Adjusted height
      )
      plot_tracks <- c(plot_tracks, list(gc_track))
    }
    
    # Plot the tracks with adjusted heights
    plotTracks(
      plot_tracks,
      from = plot_params$start,
      to = plot_params$end,
      main = paste("Genome Region for", plot_params$chromosome),
      cex.title = 1
    )
  })
}
