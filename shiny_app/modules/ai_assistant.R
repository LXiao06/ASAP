# ============================================================================
# AI Assistant Shiny Module
# ============================================================================
# Reusable module for AI chat interface that can be embedded in any tab
# ============================================================================

#' AI Assistant Module UI
#'
#' @param id Module namespace ID
#' @param height Height of the chat display area
#'
#' @return Shiny UI elements
ai_assistant_ui <- function(id, height = "400px") {
  ns <- NS(id)
  
  tagList(
    card(
      card_header(
        class = "bg-primary text-white",
        tags$i(class = "fas fa-robot me-2"),
        "AI Assistant"
      ),
      card_body(
        # Chat display area
        div(
          id = ns("chat_display"),
          style = paste0(
            "height: ", height, "; ",
            "overflow-y: auto; ",
            "border: 1px solid #dee2e6; ",
            "border-radius: 0.25rem; ",
            "padding: 1rem; ",
            "background-color: #f8f9fa; ",
            "margin-bottom: 1rem;"
          ),
          div(
            class = "text-muted text-center",
            style = "padding: 2rem;",
            tags$i(class = "fas fa-comments fa-2x mb-3"),
            tags$p("Ask me anything about ASAP functions and parameters!"),
            tags$small("Examples: 'How do I detect syllables?' or 'What threshold should I use?'")
          )
        ),
        
        # Quick action buttons
        div(
          class = "mb-2",
          actionButton(
            ns("quick_params"),
            "Suggest Parameters",
            icon = icon("sliders"),
            class = "btn-sm btn-outline-secondary me-2"
          ),
          actionButton(
            ns("quick_workflow"),
            "Workflow Help",
            icon = icon("route"),
            class = "btn-sm btn-outline-secondary me-2"
          ),
          actionButton(
            ns("clear_chat"),
            "Clear Chat",
            icon = icon("trash"),
            class = "btn-sm btn-outline-danger"
          )
        ),
        
        # Input area
        div(
          class = "input-group",
          textInput(
            ns("user_input"),
            label = NULL,
            placeholder = "Type your question here...",
            width = "100%"
          ),
          actionButton(
            ns("send_btn"),
            "Send",
            icon = icon("paper-plane"),
            class = "btn-primary"
          )
        ),
        
        # Loading indicator
        uiOutput(ns("loading_indicator"))
      )
    )
  )
}

#' AI Assistant Module Server
#'
#' @param id Module namespace ID
#' @param current_function Reactive value with current function name
#' @param current_context Reactive value with additional context
#'
#' @return Server logic for AI assistant
ai_assistant_server <- function(id, current_function = reactive(NULL), 
                                current_context = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Chat history
    chat_history <- reactiveVal(list())
    
    # Loading state
    is_loading <- reactiveVal(FALSE)
    
    # Add message to chat
    add_message <- function(role, content) {
      history <- chat_history()
      history[[length(history) + 1]] <- list(
        role = role,
        content = content,
        timestamp = Sys.time()
      )
      chat_history(history)
    }
    
    # Render chat display
    observe({
      history <- chat_history()
      
      if (length(history) == 0) {
        # Show welcome message
        chat_html <- div(
          class = "text-muted text-center",
          style = "padding: 2rem;",
          tags$i(class = "fas fa-comments fa-2x mb-3"),
          tags$p("Ask me anything about ASAP functions and parameters!"),
          tags$small("Examples: 'How do I detect syllables?' or 'What threshold should I use?'")
        )
      } else {
        # Render messages
        chat_html <- tagList(
          lapply(history, function(msg) {
            if (msg$role == "user") {
              div(
                class = "mb-3",
                div(
                  class = "d-flex justify-content-end",
                  div(
                    class = "bg-primary text-white rounded p-2",
                    style = "max-width: 80%;",
                    tags$strong(tags$i(class = "fas fa-user me-2"), "You"),
                    tags$br(),
                    msg$content
                  )
                )
              )
            } else {
              div(
                class = "mb-3",
                div(
                  class = "d-flex justify-content-start",
                  div(
                    class = "bg-white border rounded p-2",
                    style = "max-width: 80%;",
                    tags$strong(tags$i(class = "fas fa-robot me-2 text-primary"), "AI Assistant"),
                    tags$br(),
                    tags$div(
                      class = "mt-2",
                      HTML(markdown::renderMarkdown(text = msg$content))
                    )
                  )
                )
              )
            }
          })
        )
      }
      
      # Update chat display
      removeUI(paste0("#", ns("chat_display"), " > *"), multiple = TRUE)
      insertUI(
        selector = paste0("#", ns("chat_display")),
        where = "beforeEnd",
        ui = chat_html
      )
      
      # Scroll to bottom
      shinyjs::runjs(paste0(
        "document.getElementById('", ns("chat_display"), "').scrollTop = ",
        "document.getElementById('", ns("chat_display"), "').scrollHeight;"
      ))
    })
    
    # Loading indicator
    output$loading_indicator <- renderUI({
      if (is_loading()) {
        div(
          class = "text-center mt-2",
          tags$i(class = "fas fa-spinner fa-spin me-2"),
          tags$small("AI is thinking...")
        )
      }
    })
    
    # Send message
    send_message <- function(message) {
      if (message == "" || is.null(message)) return()
      
      # Add user message
      add_message("user", message)
      
      # Clear input
      updateTextInput(session, "user_input", value = "")
      
      # Set loading state
      is_loading(TRUE)
      
      # Get AI response
      future::future({
        # Build context
        context_parts <- c()
        if (!is.null(current_function())) {
          context_parts <- c(context_parts, 
                           paste("Current function:", current_function()))
        }
        if (!is.null(current_context())) {
          context_parts <- c(context_parts, current_context())
        }
        context <- if (length(context_parts) > 0) {
          paste(context_parts, collapse = "\n")
        } else {
          NULL
        }
        
        # Call AI
        get_ai_suggestion(message, context = context)
      }) %...>% {
        # Add AI response
        add_message("assistant", .)
        is_loading(FALSE)
      }
    }
    
    # Send button click
    observeEvent(input$send_btn, {
      send_message(input$user_input)
    })
    
    # Enter key press
    observeEvent(input$user_input, {
      if (grepl("\n$", input$user_input)) {
        send_message(trimws(input$user_input))
      }
    })
    
    # Quick action: Suggest parameters
    observeEvent(input$quick_params, {
      func <- current_function()
      if (is.null(func)) {
        send_message("What function would you like parameter suggestions for?")
      } else {
        send_message(paste("What parameters should I use for", func, "?"))
      }
    })
    
    # Quick action: Workflow help
    observeEvent(input$quick_workflow, {
      send_message("Can you explain the typical workflow for motif detection in ASAP?")
    })
    
    # Clear chat
    observeEvent(input$clear_chat, {
      chat_history(list())
      clear_ai_cache()
      show_notification("Chat cleared", "info")
    })
    
    # Return reactive values
    return(list(
      chat_history = chat_history
    ))
  })
}
