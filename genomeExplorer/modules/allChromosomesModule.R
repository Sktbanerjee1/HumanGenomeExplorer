# Module for rendering the "All Human Chromosomes" plot
allChromosomesModule <- function(input, output, session) {
  
  output$allChromosomesPlot <- renderPlot({
    
    # If the RDS file does not exist, generate the plot
    print("Generating plot...")
    chromosome_list <- paste0("chr", c(1:22, "X", "Y"))
    
    # Create ideogram tracks for all chromosomes and explicitly set the plotting range
    plot_list <- lapply(chromosome_list, function(chrom) {
      ideogram_track <- IdeogramTrack(
        genome = "hg38",
        chromosome = chrom,
        showBandId = FALSE,
        showId = FALSE,
        useGenomeCoords = TRUE
      )
      
      from <- 1
      to <- 1
      
      # Capture the plot as a grob with explicit 'from' and 'to' arguments
      grid.grabExpr(plotTracks(ideogram_track, from = from, to = to, main = chrom))
    })
    
    # Arrange the plots in a grid layout
    final_plot <- do.call(gridExtra::grid.arrange, c(plot_list, ncol = 3))
    
    # Return the final plot
    final_plot
  })
}
