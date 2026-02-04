# ============================================================================
# Main Server Logic for ASAP Shiny App
# ============================================================================

server_main <- function(input, output, session) {
  
  # Initialize AI assistant for home page
  ai_assistant_server("home_ai")
  
  # Reactive values for app state
  app_state <- reactiveValues(
    current_wav = NULL,
    original_filename = NULL,  # Store original filename
    wav_info = NULL,
    bout_results = NULL,
    syllable_results = NULL,
    bout_plot = NULL,
    segment_plot = NULL,
    last_output = "none"  # Track which output to show: "viz", "bout", "segment"
  )
  
  # Navigate to Single WAV Analysis when "Get Started" is clicked
  observeEvent(input$get_started, {
    updateNavbarPage(session, inputId = "navbar", selected = "single_wav")
  })
  
  # Load Single WAV Analysis UI
  output$single_wav_content <- renderUI({
    source("ui/ui_single_wav.R", local = TRUE)
    ui_single_wav()
  })
  
  # Load Motif Detection UI
  output$motif_detection_content <- renderUI({
    card(
      card_header("Motif Detection"),
      card_body(
        h4("Coming Soon"),
        p("The motif detection workflow is under development."),
        p("This will include:"),
        tags$ul(
          tags$li("Interactive template creation"),
          tags$li("Template matching and detection"),
          tags$li("Motif extraction and visualization"),
          tags$li("AI-guided parameter optimization")
        )
      )
    )
  })
  
  # Load Longitudinal Analysis UI
  output$longitudinal_content <- renderUI({
    card(
      card_header("Longitudinal Analysis"),
      card_body(
        h4("Coming Soon"),
        p("The longitudinal analysis workflow is under development."),
        p("This will include:"),
        tags$ul(
          tags$li("SAP object creation"),
          tags$li("Bulk processing pipelines"),
          tags$li("Heatmap visualizations"),
          tags$li("UMAP and clustering analysis")
        )
      )
    )
  })
  
  # Single WAV Analysis Server Logic
  # ============================================================================
  
  # Initialize AI assistant for single WAV
  single_wav_ai <- ai_assistant_server(
    "single_wav_ai",
    current_function = reactive({
      # Return current function based on active analysis
      if (!is.null(input$run_bout_detection)) "find_bout"
      else if (!is.null(input$run_segmentation)) "segment"
      else NULL
    })
  )
  
  # File upload handling
  observeEvent(input$wav_file, {
    req(input$wav_file)
    
    file_path <- input$wav_file$datapath
    validation <- validate_wav_file(file_path)
    
    if (validation$valid) {
      app_state$current_wav <- file_path
      app_state$original_filename <- input$wav_file$name  # Store original filename
      app_state$wav_info <- validation
      
      # Update time range sliders based on file duration
      duration <- validation$duration
      
      # Update visualization time range
      updateSliderInput(
        session,
        "viz_time_range",
        max = duration,
        value = c(0, min(10, duration))
      )
      
      # Update segmentation time range slider
      updateSliderInput(
        session,
        "seg_time_range",
        max = duration,
        value = c(min(1, duration), min(5, duration))
      )
      
      show_notification(
        paste("File loaded successfully:", input$wav_file$name),
        "success"
      )
    } else {
      show_notification(validation$message, "error")
    }
  })
  
  # Display file info
  output$file_info <- renderUI({
    req(app_state$wav_info)
    info <- app_state$wav_info
    
    div(
      class = "alert alert-info p-2",
      style = "font-size: 0.85rem;",
      tags$strong("File Info:"),
      tags$ul(
        class = "mb-0 mt-1",
        style = "padding-left: 1.2rem;",
        tags$li(paste("Duration:", sprintf("%.2f", info$duration), "s")),
        tags$li(paste("Sample Rate:", info$sample_rate, "Hz")),
        tags$li(paste("Channels:", info$channels))
      )
    )
  })
  
  # Output title (dynamic based on last action)
  output$output_title <- renderText({
    switch(app_state$last_output,
      "viz" = "Spectrogram Visualization",
      "bout" = "Bout Detection Results",
      "segment" = "Syllable Segmentation Results",
      "Output Window"
    )
  })
  
  # Unified main output plot
  output$main_output_plot <- renderPlot({
    if (app_state$last_output == "none") {
      # Show welcome message
      plot.new()
      text(0.5, 0.5, "Upload a WAV file and run an analysis to see results here", 
           cex = 1.2, col = "gray50")
      return()
    }
    
    if (app_state$last_output == "viz") {
      req(app_state$current_wav)
      tryCatch({
        visualize_song(
          app_state$current_wav,
          start_time_in_second = input$viz_time_range[1],
          end_time_in_second = input$viz_time_range[2]
        )
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error:", conditionMessage(e)), col = "red")
      })
    } else if (app_state$last_output == "bout") {
      req(app_state$bout_plot)
      print(app_state$bout_plot)
    } else if (app_state$last_output == "segment") {
      req(app_state$segment_plot)
      print(app_state$segment_plot)
    }
  })
  
  # Visualization button
  observeEvent(input$run_visualization, {
    req(app_state$current_wav)
    app_state$last_output <- "viz"
  })
  
  # View full recording
  observeEvent(input$viz_full, {
    req(app_state$wav_info)
    updateSliderInput(
      session,
      "viz_time_range",
      value = c(0, app_state$wav_info$duration)
    )
    app_state$last_output <- "viz"
  })
  
  # Bout detection
  observeEvent(input$run_bout_detection, {
    req(app_state$current_wav)
    
    withProgress(message = "Detecting bouts...", {
      tryCatch({
        bouts <- find_bout(
          app_state$current_wav,
          rms_threshold = input$bout_rms_threshold,
          min_duration = input$bout_min_duration,
          plot = TRUE
        )
        
        app_state$bout_results <- bouts
        
        # Capture the plot
        app_state$bout_plot <- recordPlot()
        app_state$last_output <- "bout"
        
        if (!is.null(bouts) && nrow(bouts) > 0) {
          show_notification(
            paste("Found", nrow(bouts), "bout(s)"),
            "success"
          )
        } else {
          show_notification(
            "No bouts detected. Try adjusting parameters.",
            "warning"
          )
        }
      }, error = function(e) {
        show_notification(
          paste("Error:", conditionMessage(e)),
          "error"
        )
      })
    })
  })
  
  # Display bout results
  output$bout_results <- renderUI({
    req(app_state$bout_results)
    
    div(
      class = "alert alert-success p-2",
      style = "font-size: 0.85rem;",
      tags$strong(paste("Found", nrow(app_state$bout_results), "bout(s)")),
      tags$br(),
      tags$small("See output window above")
    )
  })
  
  # Syllable segmentation
  observeEvent(input$run_segmentation, {
    req(app_state$current_wav)
    
    withProgress(message = "Segmenting syllables...", {
      tryCatch({
        # Run segmentation with plot=TRUE to use original visualization
        syllables <- segment(
          app_state$current_wav,
          start_time = input$seg_time_range[1],
          end_time = input$seg_time_range[2],
          flim = input$seg_flim,
          silence_threshold = input$seg_silence_threshold,
          min_syllable_ms = input$seg_syllable_duration_range[1],
          max_syllable_ms = input$seg_syllable_duration_range[2],
          min_level_db = input$seg_level_range[1],
          max_level_db = input$seg_level_range[2],
          verbose = FALSE,
          plot = TRUE  # Use original plot with boxes
        )
        
        # Capture the plot that was just created
        app_state$segment_plot <- recordPlot()
        
        app_state$syllable_results <- syllables
        app_state$last_output <- "segment"
        
        if (!is.null(syllables) && nrow(syllables) > 0) {
          show_notification(
            paste("Found", nrow(syllables), "syllable(s)"),
            "success"
          )
        } else {
          show_notification(
            "No syllables detected. Try adjusting parameters.",
            "warning"
          )
        }
      }, error = function(e) {
        show_notification(
          paste("Error:", conditionMessage(e)),
          "error"
        )
      })
    })
  })
  
  # Display segmentation results
  output$segmentation_results <- renderUI({
    req(app_state$syllable_results)
    
    div(
      class = "alert alert-success p-2",
      style = "font-size: 0.85rem;",
      tags$strong(paste("Found", nrow(app_state$syllable_results), "syllable(s)")),
      tags$br(),
      tags$small("See output window above")
    )
  })
  
  # Results table (shows latest results)
  output$results_table <- renderDT({
    if (!is.null(app_state$syllable_results) && app_state$last_output == "segment") {
      # Format syllable results
      results <- app_state$syllable_results
      
      # Replace filename with original
      if ("filename" %in% names(results) && !is.null(app_state$original_filename)) {
        results$filename <- app_state$original_filename
      }
      
      # Format numeric columns to 2 decimals
      numeric_cols <- sapply(results, is.numeric)
      results[numeric_cols] <- lapply(results[numeric_cols], function(x) round(x, 2))
      
      datatable(
        results,
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE,
        caption = "Syllable Segmentation Results"
      )
    } else if (!is.null(app_state$bout_results) && app_state$last_output == "bout") {
      # Format bout results
      results <- app_state$bout_results
      
      # Replace filename with original
      if ("filename" %in% names(results) && !is.null(app_state$original_filename)) {
        results$filename <- app_state$original_filename
      }
      
      # Format numeric columns to 2 decimals
      numeric_cols <- sapply(results, is.numeric)
      results[numeric_cols] <- lapply(results[numeric_cols], function(x) round(x, 2))
      
      datatable(
        results,
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE,
        caption = "Bout Detection Results"
      )
    } else {
      datatable(
        data.frame(Message = "No results yet. Run an analysis to see results."),
        options = list(dom = 't'),
        rownames = FALSE
      )
    }
  })
}
