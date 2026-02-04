# ============================================================================
# Single WAV Analysis UI - Redesigned Layout
# ============================================================================

ui_single_wav <- function() {
  tagList(
    # Main layout: Left sidebar (functions) + Right output area
    layout_sidebar(
      sidebar = sidebar(
        width = 450,  # Wider sidebar
        open = "always",
        
        # File upload
        card(
          card_header(
            icon("upload"),
            "Upload Audio File"
          ),
          card_body(
            class = "p-2",
            fileInput(
              "wav_file",
              NULL,
              accept = c(".wav", ".WAV"),
              buttonLabel = "Browse...",
              placeholder = "No file selected"
            ),
            uiOutput("file_info")
          )
        ),
        
        # Visualization controls
        card(
          card_header(
            icon("chart-area"),
            "Visualization"
          ),
          card_body(
            class = "p-2",
            sliderInput(
              "viz_time_range",
              "Time Range (seconds)",
              min = 0,
              max = 10,
              value = c(0, 10),
              step = 0.1,
              width = "100%",
              ticks = FALSE  # Remove secondary scale
            ),
            actionButton(
              "run_visualization",
              "View Spectrogram",
              icon = icon("play"),
              class = "btn-primary w-100 mb-2"
            ),
            actionButton(
              "viz_full",
              "View Full Recording",
              icon = icon("expand"),
              class = "btn-sm btn-outline-secondary w-100"
            )
          )
        ),
        
        # Bout detection
        card(
          card_header(
            icon("wave-square"),
            "Bout Detection",
            popover(
              icon("circle-question"),
              "Detect continuous singing periods"
            )
          ),
          card_body(
            class = "p-2",
            sliderInput(
              "bout_rms_threshold",
              "RMS Threshold",
              min = 0,
              max = 1,
              value = DEFAULT_PARAMS$rms_threshold,
              step = 0.01,
              ticks = FALSE  # Remove secondary scale
            ),
            sliderInput(
              "bout_min_duration",
              "Min Duration (s)",
              min = 0.1,
              max = 5,
              value = DEFAULT_PARAMS$min_duration,
              step = 0.1,
              ticks = FALSE
            ),
            actionButton(
              "run_bout_detection",
              "Detect Bouts",
              icon = icon("play"),
              class = "btn-primary w-100"
            ),
            hr(),
            uiOutput("bout_results")
          )
        ),
        
        # Syllable segmentation
        card(
          card_header(
            icon("scissors"),
            "Syllable Segmentation",
            popover(
              icon("circle-question"),
              "Detect individual syllables"
            )
          ),
          card_body(
            class = "p-2",
            sliderInput(
              "seg_time_range",
              "Time Range (s)",
              min = 0,
              max = 10,
              value = c(1, 5),
              step = 0.1,
              ticks = FALSE
            ),
            sliderInput(
              "seg_flim",
              "Frequency Limits (kHz)",
              min = 0,
              max = 20,
              value = c(1, 8),
              step = 0.5,
              ticks = FALSE
            ),
            sliderInput(
              "seg_silence_threshold",
              "Silence Threshold",
              min = 0,
              max = 0.1,
              value = DEFAULT_PARAMS$silence_threshold,
              step = 0.001,
              ticks = FALSE
            ),
            sliderInput(
              "seg_syllable_duration_range",
              "Syllable Duration (ms)",
              min = 10,
              max = 500,
              value = c(DEFAULT_PARAMS$min_syllable_ms, DEFAULT_PARAMS$max_syllable_ms),
              step = 5,
              ticks = FALSE
            ),
            sliderInput(
              "seg_level_range",
              "Level Range (dB)",
              min = 0,
              max = 50,
              value = c(10, 30),
              step = 1,
              ticks = FALSE
            ),
            actionButton(
              "run_segmentation",
              "Segment Syllables",
              icon = icon("play"),
              class = "btn-primary w-100"
            ),
            hr(),
            uiOutput("segmentation_results")
          )
        )
      ),
      
      # Main output area (right side)
      layout_column_wrap(
        width = 1,
        heights_equal = "row",
        
        # Unified output window - TALLER
        card(
          full_screen = TRUE,
          height = "700px",  # Increased from 600px
          card_header(
            icon("chart-bar"),
            textOutput("output_title", inline = TRUE)
          ),
          card_body(
            plotOutput("main_output_plot", height = "620px")  # Increased from 500px
          )
        ),
        
        # Results table
        card(
          card_header(
            icon("table"),
            "Analysis Results"
          ),
          card_body(
            DTOutput("results_table")
          )
        ),
        
        # AI Assistant at bottom
        ai_assistant_ui("single_wav_ai", height = "400px")
      )
    )
  )
}
