# ============================================================================
# Test Script for ASAP Shiny App
# ============================================================================
# This script checks dependencies and tests core functionality before
# launching the full app.
# ============================================================================

cat("Testing ASAP Shiny App Setup...\n\n")

# Check required packages -------------------------------------------------
cat("1. Checking required packages...\n")

required_pkgs <- c(
  "shiny", "bslib", "ASAP", "httr2", "jsonlite",
  "plotly", "DT", "shinyWidgets", "digest", "markdown",
  "tuneR", "monitoR"
)

missing_pkgs <- c()
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, pkg)
    cat("   ✗", pkg, "- NOT INSTALLED\n")
  } else {
    cat("   ✓", pkg, "\n")
  }
}

if (length(missing_pkgs) > 0) {
  cat("\n❌ Missing packages detected!\n")
  cat("Install them with:\n")
  cat("install.packages(c(", paste0("'", missing_pkgs, "'", collapse = ", "), "))\n\n")
  stop("Please install missing packages before running the app.")
} else {
  cat("\n✅ All required packages are installed!\n\n")
}

# Check API key configuration ---------------------------------------------
cat("2. Checking API key configuration...\n")

api_key <- Sys.getenv("OPENAI_API_KEY")
if (api_key == "" || api_key == "your-api-key-here") {
  cat("   ⚠️  OpenAI API key not configured\n")
  cat("   The app will run but AI features will be disabled.\n")
  cat("   To enable AI features:\n")
  cat("   1. Copy .Renviron.example to .Renviron\n")
  cat("   2. Add your OpenAI API key\n")
  cat("   3. Restart R session\n\n")
} else {
  cat("   ✓ API key configured\n")
  cat("   Model:", Sys.getenv("OPENAI_MODEL", "gpt-3.5-turbo"), "\n\n")
}

# Check example data ------------------------------------------------------
cat("3. Checking example data...\n")

example_wav <- system.file("extdata", "zf_example.wav", package = "ASAP")
if (file.exists(example_wav)) {
  cat("   ✓ Example WAV file found\n")
  cat("   Path:", example_wav, "\n\n")
} else {
  cat("   ⚠️  Example WAV file not found\n")
  cat("   You can still use the app with your own WAV files.\n\n")
}

# Test AI integration (if configured) -------------------------------------
if (api_key != "" && api_key != "your-api-key-here") {
  cat("4. Testing AI integration...\n")
  
  tryCatch({
    source("shiny_app/utils/ai_integration.R")
    
    cat("   Testing API connection...\n")
    response <- get_ai_suggestion("Say 'Hello' in one word.", use_cache = FALSE)
    
    if (grepl("Error|❌", response)) {
      cat("   ✗ AI test failed\n")
      cat("   Response:", response, "\n\n")
    } else {
      cat("   ✓ AI integration working!\n")
      cat("   Response:", response, "\n\n")
    }
  }, error = function(e) {
    cat("   ✗ Error testing AI:", conditionMessage(e), "\n\n")
  })
} else {
  cat("4. Skipping AI integration test (no API key)\n\n")
}

# Summary -----------------------------------------------------------------
cat("=" , rep("=", 70), "\n", sep = "")
cat("SETUP SUMMARY\n")
cat("=" , rep("=", 70), "\n", sep = "")

if (length(missing_pkgs) == 0) {
  cat("✅ All dependencies installed\n")
} else {
  cat("❌ Missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
}

if (api_key != "" && api_key != "your-api-key-here") {
  cat("✅ AI features enabled\n")
} else {
  cat("⚠️  AI features disabled (no API key)\n")
}

if (file.exists(example_wav)) {
  cat("✅ Example data available\n")
} else {
  cat("⚠️  No example data (use your own WAV files)\n")
}

cat("\n")

if (length(missing_pkgs) == 0) {
  cat("🚀 Ready to launch!\n")
  cat("Run the app with: shiny::runApp('shiny_app')\n\n")
} else {
  cat("⚠️  Please install missing packages first.\n\n")
}
