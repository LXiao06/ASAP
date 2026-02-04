# ============================================================================
# Main UI for ASAP Shiny App
# ============================================================================

# Source module files
source("modules/ai_assistant.R", local = TRUE)

# Main UI function
ui_main <- function() {
  page_navbar(
    title = div(
      tags$img(
        src = "logo.png",
        height = "30px",
        class = "me-2",
        onerror = "this.style.display='none'"
      ),
      "ASAP Studio"
    ),
    theme = app_theme,
    fillable = TRUE,
    
    # Home / Welcome Tab
    nav_panel(
      title = "Home",
      icon = icon("home"),
      layout_columns(
        col_widths = c(8, 4),
        
        # Main content
        card(
          card_header("Welcome to ASAP Studio"),
          card_body(
            h4("Automated Sound Analysis Pipeline"),
            p(
              "This interactive application provides a user-friendly interface for",
              "analyzing zebra finch vocalizations and other bird songs using the ASAP R package."
            ),
            
            h5("Features", class = "mt-4"),
            tags$ul(
              tags$li(
                tags$strong("AI-Powered Assistant:"),
                "Get intelligent suggestions for parameters and workflows"
              ),
              tags$li(
                tags$strong("Interactive Visualizations:"),
                "Explore spectrograms and analysis results"
              ),
              tags$li(
                tags$strong("Guided Workflows:"),
                "Step-by-step pipelines for common analysis tasks"
              ),
              tags$li(
                tags$strong("Real-time Feedback:"),
                "Immediate validation and optimization suggestions"
              )
            ),
            
            h5("Quick Start", class = "mt-4"),
            tags$ol(
              tags$li("Upload a WAV file in the 'Single WAV Analysis' tab"),
              tags$li("Visualize your audio with spectrograms"),
              tags$li("Use the AI assistant for parameter guidance"),
              tags$li("Run analysis and explore results")
            ),
            
            div(
              class = "mt-4",
              actionButton(
                "get_started",
                "Get Started with Single WAV Analysis",
                icon = icon("play"),
                class = "btn-primary btn-lg"
              )
            )
          )
        ),
        
        # AI Assistant sidebar
        ai_assistant_ui("home_ai", height = "500px")
      )
    ),
    
    # Single WAV Analysis Tab
    nav_panel(
      title = "Single WAV Analysis",
      icon = icon("file-audio"),
      value = "single_wav",
      # Content will be loaded from ui_single_wav.R
      uiOutput("single_wav_content")
    ),
    
    # Motif Detection Tab
    nav_panel(
      title = "Motif Detection",
      icon = icon("search"),
      value = "motif_detection",
      # Content will be loaded from ui_motif_detection.R
      uiOutput("motif_detection_content")
    ),
    
    # Longitudinal Analysis Tab
    nav_panel(
      title = "Longitudinal Analysis",
      icon = icon("chart-line"),
      value = "longitudinal",
      # Content will be loaded from ui_longitudinal.R
      uiOutput("longitudinal_content")
    ),
    
    # Settings Tab
    nav_panel(
      title = "Settings",
      icon = icon("cog"),
      card(
        card_header("Application Settings"),
        card_body(
          h5("AI Configuration"),
          textInput(
            "ai_model",
            "OpenAI Model",
            value = AI_CONFIG$model,
            placeholder = "gpt-3.5-turbo"
          ),
          sliderInput(
            "ai_temperature",
            "AI Temperature (creativity)",
            min = 0,
            max = 1,
            value = AI_CONFIG$temperature,
            step = 0.1
          ),
          
          hr(),
          
          h5("Default Parameters"),
          p("Set default values for ASAP analysis functions"),
          
          numericInput(
            "default_freq_min",
            "Default Frequency Min (kHz)",
            value = DEFAULT_PARAMS$template_freq_min,
            min = 0,
            max = 20,
            step = 0.5
          ),
          numericInput(
            "default_freq_max",
            "Default Frequency Max (kHz)",
            value = DEFAULT_PARAMS$template_freq_max,
            min = 0,
            max = 20,
            step = 0.5
          ),
          
          hr(),
          
          h5("About"),
          p(
            "ASAP Studio v1.0.0",
            tags$br(),
            "Automated Sound Analysis Pipeline",
            tags$br(),
            "Built with Shiny and bslib",
            tags$br(),
            tags$a(
              href = "https://github.com/LXiao06/ASAP",
              target = "_blank",
              "View on GitHub"
            )
          )
        )
      )
    ),
    
    # Footer
    nav_spacer(),
    nav_item(
      tags$a(
        icon("github"),
        "GitHub",
        href = "https://github.com/LXiao06/ASAP",
        target = "_blank",
        class = "nav-link"
      )
    )
  )
}
