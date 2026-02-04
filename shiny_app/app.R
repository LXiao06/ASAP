# ============================================================================
# ASAP Studio - Main Application File
# ============================================================================
# This is the main entry point for the ASAP Studio interactive analysis app.
# It loads the UI and server components and launches the app.
#
# To run this app:
# 1. Set up your .Renviron file with OPENAI_API_KEY
# 2. Install required packages (see README.md)
# 3. Run: shiny::runApp("shiny_app")
# ============================================================================

# Load global configuration
source("global.R", local = TRUE)

# Load UI components
source("ui/ui_main.R", local = TRUE)

# Load server components
source("server/server_main.R", local = TRUE)

# Create the Shiny app
shinyApp(
  ui = ui_main(),
  server = server_main
)
