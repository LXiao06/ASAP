# ============================================================================
# AI Integration Utilities for ASAP Shiny App
# ============================================================================
# This module provides functions for integrating OpenAI API to assist users
# with parameter selection and workflow guidance in the ASAP package.
#
# Author: ASAP Development Team
# Date: 2026-02-03
# ============================================================================

library(httr2)
library(jsonlite)

# Global AI configuration
AI_CONFIG <- list(
  api_base = "https://api.openai.com/v1",
  model = Sys.getenv("OPENAI_MODEL", "gpt-3.5-turbo"),
  max_tokens = 500,
  temperature = 0.7,
  timeout = 30
)

# Cache for AI responses to reduce API calls
ai_response_cache <- new.env()

#' Get AI Suggestion
#'
#' Sends a user question to OpenAI API and returns an AI-generated response
#' with context about the ASAP package.
#'
#' @param user_question Character string with the user's question
#' @param context Optional additional context to include in the prompt
#' @param use_cache Logical, whether to use cached responses (default: TRUE)
#'
#' @return Character string with AI response, or error message if API call fails
#'
#' @examples
#' \dontrun{
#' response <- get_ai_suggestion(
#'   "What parameters should I use for detect_template?"
#' )
#' }
get_ai_suggestion <- function(user_question, context = NULL, use_cache = TRUE) {
  
  # Check for API key
  api_key <- Sys.getenv("OPENAI_API_KEY")
  if (api_key == "" || api_key == "your-api-key-here") {
    return(paste(
      "⚠️ OpenAI API key not configured.",
      "Please set OPENAI_API_KEY in your .Renviron file.",
      "\nSee .Renviron.example for instructions."
    ))
  }
  
  # Check cache
  cache_key <- digest::digest(paste(user_question, context))
  if (use_cache && exists(cache_key, envir = ai_response_cache)) {
    return(get(cache_key, envir = ai_response_cache))
  }
  
  # Build system prompt with ASAP context
  system_prompt <- paste(
    "You are an expert assistant for the ASAP (Acoustic Signal Analysis Pipeline) R package,",
    "which is used for analyzing zebra finch vocalizations and other bird songs.",
    "Your role is to help users understand functions, choose appropriate parameters,",
    "and optimize their audio analysis workflows.",
    "\n\nKey ASAP workflows:",
    "1. Single WAV analysis: visualize_song(), find_bout(), segment()",
    "2. Motif detection: create_audio_clip(), create_template(), detect_template(), find_motif()",
    "3. Longitudinal analysis: create_sap_object(), analyze_spectral(), find_clusters(), run_umap()",
    "\n\nProvide concise, practical advice. Focus on parameter values and workflow steps.",
    "Use bullet points for clarity. Mention typical values for zebra finch analysis when relevant."
  )
  
  # Add additional context if provided
  if (!is.null(context)) {
    system_prompt <- paste(system_prompt, "\n\nAdditional context:", context)
  }
  
  # Build messages
  messages <- list(
    list(role = "system", content = system_prompt),
    list(role = "user", content = user_question)
  )
  
  # Make API request with error handling
  tryCatch({
    response <- request(paste0(AI_CONFIG$api_base, "/chat/completions")) |>
      req_headers(
        "Authorization" = paste("Bearer", api_key),
        "Content-Type" = "application/json"
      ) |>
      req_body_json(list(
        model = AI_CONFIG$model,
        messages = messages,
        temperature = AI_CONFIG$temperature,
        max_tokens = AI_CONFIG$max_tokens
      )) |>
      req_timeout(AI_CONFIG$timeout) |>
      req_retry(max_tries = 3) |>
      req_perform() |>
      resp_body_json()
    
    # Extract response text
    ai_response <- response$choices[[1]]$message$content
    
    # Cache the response
    assign(cache_key, ai_response, envir = ai_response_cache)
    
    return(ai_response)
    
  }, error = function(e) {
    error_msg <- paste(
      "❌ Error calling OpenAI API:",
      conditionMessage(e),
      "\n\nPlease check:",
      "- Your API key is valid",
      "- You have sufficient API credits",
      "- Your internet connection is working"
    )
    return(error_msg)
  })
}

#' Get Parameter Help
#'
#' Specialized function to get AI assistance for function parameters.
#'
#' @param function_name Name of the ASAP function
#' @param user_context Description of what the user is trying to do
#'
#' @return Character string with parameter recommendations
#'
#' @examples
#' \dontrun{
#' help <- get_parameter_help(
#'   function_name = "detect_template",
#'   user_context = "detecting syllables in zebra finch songs"
#' )
#' }
get_parameter_help <- function(function_name, user_context = "") {
  
  # Build specific question about parameters
  question <- paste0(
    "I want to use the `", function_name, "()` function in ASAP",
    if (user_context != "") paste(" for", user_context) else "",
    ".\n\nPlease tell me:",
    "\n1. Which parameters are REQUIRED?",
    "\n2. Which parameters are optional but ESSENTIAL for good results?",
    "\n3. What are typical values for zebra finch analysis?",
    "\n\nFormat your response with clear sections and bullet points."
  )
  
  # Get function documentation as context
  context <- tryCatch({
    # Try to get function help
    help_text <- capture.output(help(function_name, package = "ASAP"))
    if (length(help_text) > 0) {
      paste("Function documentation available:", 
            paste(head(help_text, 20), collapse = "\n"))
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  return(get_ai_suggestion(question, context = context))
}

#' Validate Parameters with AI
#'
#' Use AI to validate parameter combinations and suggest improvements.
#'
#' @param function_name Name of the ASAP function
#' @param params Named list of parameters and their values
#'
#' @return Character string with validation feedback
#'
#' @examples
#' \dontrun{
#' feedback <- validate_parameters(
#'   function_name = "segment",
#'   params = list(
#'     silence_threshold = 0.01,
#'     min_syllable_ms = 20,
#'     max_syllable_ms = 240
#'   )
#' )
#' }
validate_parameters <- function(function_name, params) {
  
  # Format parameters for the question
  param_text <- paste(
    sapply(names(params), function(p) {
      paste0("- ", p, " = ", params[[p]])
    }),
    collapse = "\n"
  )
  
  question <- paste0(
    "I'm using `", function_name, "()` with these parameters:\n\n",
    param_text,
    "\n\nPlease review these parameters and:",
    "\n1. Identify any potential issues or conflicts",
    "\n2. Suggest improvements if needed",
    "\n3. Confirm if this looks good for zebra finch analysis",
    "\n\nBe concise and specific."
  )
  
  return(get_ai_suggestion(question, use_cache = FALSE))
}

#' Build ASAP Context for AI
#'
#' Retrieves function documentation and examples to provide context to AI.
#'
#' @param function_name Name of the ASAP function
#'
#' @return Character string with formatted context
build_asap_context <- function(function_name) {
  
  context_parts <- c()
  
  # Try to get function documentation
  tryCatch({
    help_file <- utils:::.getHelpFile(help(function_name, package = "ASAP"))
    if (!is.null(help_file)) {
      context_parts <- c(context_parts, 
                        "Function documentation available from ASAP package")
    }
  }, error = function(e) {})
  
  # Add common parameter ranges for known functions
  common_params <- list(
    detect_template = "threshold: 0.3-0.7, proximity_window: 0.5-2.0 seconds",
    create_template = "freq_min: 1 kHz, freq_max: 8-10 kHz for zebra finch",
    segment = "min_syllable_ms: 20-50, max_syllable_ms: 200-300",
    find_bout = "rms_threshold: 0.05-0.15, min_duration: 0.5-1.0 seconds"
  )
  
  if (function_name %in% names(common_params)) {
    context_parts <- c(context_parts, 
                      paste("Typical parameters:", common_params[[function_name]]))
  }
  
  return(paste(context_parts, collapse = "\n"))
}

#' Clear AI Response Cache
#'
#' Clears the cached AI responses. Useful for testing or when you want
#' fresh responses.
clear_ai_cache <- function() {
  rm(list = ls(envir = ai_response_cache), envir = ai_response_cache)
  message("AI response cache cleared")
}

#' Get AI Token Usage Estimate
#'
#' Estimates token usage for a given text (rough approximation).
#'
#' @param text Character string to estimate tokens for
#'
#' @return Approximate token count
estimate_tokens <- function(text) {
  # Rough estimate: ~4 characters per token
  ceiling(nchar(text) / 4)
}
