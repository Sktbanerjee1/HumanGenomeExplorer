library(shiny)
library(shinydashboard)
library(shinydashboardPlus)

# Load the UI and server functions from separate files
source("ui.R")
source("server.R")

# Run the Shiny app
shinyApp(ui = ui, server = server)
