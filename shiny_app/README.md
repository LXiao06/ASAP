# ASAP Studio

Interactive web application for the ASAP (Automated Sound Analysis Pipeline) R package with AI-powered assistance.

## Features

- 🎵 **Single WAV Analysis**: Upload and analyze individual audio files
- 🔍 **Motif Detection**: Template-based motif detection with interactive parameter optimization
- 📊 **Longitudinal Analysis**: Bulk processing and visualization of longitudinal recordings
- 🤖 **AI Assistant**: Intelligent parameter suggestions and workflow guidance
- 📈 **Interactive Visualizations**: Spectrograms, heatmaps, and UMAP plots

## Installation

### Prerequisites

1. **R** (version ≥ 4.0.0)
2. **ASAP Package** installed
3. **OpenAI API Key** (for AI assistant features)

### Required R Packages

Install the required packages:

```r
install.packages(c(
  "shiny",
  "bslib",
  "httr2",
  "jsonlite",
  "plotly",
  "DT",
  "shinyWidgets",
  "digest",
  "markdown",
  "future",
  "promises"
))
```

## Setup

### 1. Configure API Key

Copy the example environment file:

```bash
cp .Renviron.example .Renviron
```

Edit `.Renviron` and add your OpenAI API key:

```
OPENAI_API_KEY=sk-proj-your-actual-api-key-here
OPENAI_MODEL=gpt-3.5-turbo
```

**Get an API key:**
1. Go to https://platform.openai.com/api-keys
2. Sign up or log in
3. Create a new API key
4. Copy and paste it into your `.Renviron` file

**Important:** Never commit your `.Renviron` file to git! It's already in `.gitignore`.

### 2. Verify Installation

Test that all packages are installed:

```r
# In R console
required_pkgs <- c("shiny", "bslib", "ASAP", "httr2", "jsonlite", 
                   "plotly", "DT", "shinyWidgets", "digest")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  cat("Missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
} else {
  cat("All required packages are installed!\n")
}
```

## Running the App

### From RStudio

1. Open `app.R` in RStudio
2. Click the "Run App" button
3. The app will open in your browser or RStudio viewer

### From R Console

```r
# Navigate to the ASAP repository
setwd("/path/to/ASAP")

# Run the app
shiny::runApp("shiny_app")
```

### Specify Port

```r
shiny::runApp("shiny_app", port = 8080)
```

## Usage Guide

### Single WAV Analysis

1. **Upload File**: Click "Browse" and select a WAV file
2. **Visualize**: View the spectrogram and adjust time range
3. **Detect Bouts**: Set parameters and click "Detect Bouts"
4. **Segment Syllables**: Define time window and run segmentation
5. **Ask AI**: Use the AI assistant for parameter suggestions

### AI Assistant

The AI assistant can help with:
- Parameter selection for any ASAP function
- Workflow guidance
- Troubleshooting analysis issues
- Explaining function behavior

**Example questions:**
- "What parameters should I use for detect_template?"
- "How do I optimize syllable segmentation?"
- "What's a typical RMS threshold for zebra finch?"

### Tips

- Start with default parameters and adjust based on results
- Use the AI assistant's "Suggest Parameters" button for quick help
- Visualize your data before running analysis
- Check the "Quick Tips" panel for workflow guidance

## Project Structure

```
shiny_app/
├── app.R                 # Main app entry point
├── global.R              # Global configuration
├── ui/                   # UI components
│   ├── ui_main.R
│   └── ui_single_wav.R
├── server/               # Server logic
│   └── server_main.R
├── modules/              # Reusable modules
│   └── ai_assistant.R
├── utils/                # Helper functions
│   └── ai_integration.R
├── www/                  # Static assets
│   └── custom.css
├── .Renviron.example     # API key template
└── README.md             # This file
```

## Troubleshooting

### "API key not configured" error

Make sure you've created a `.Renviron` file with your OpenAI API key:

```
OPENAI_API_KEY=your-actual-key-here
```

Then restart your R session.

### "Package not found" error

Install missing packages:

```r
install.packages("package_name")
```

### WAV file won't upload

- Check that the file is a valid WAV format
- Ensure file size is under 100 MB
- Try converting the file with audio software if needed

### Slow performance

- Use `gpt-3.5-turbo` instead of `gpt-4` for faster AI responses
- Process smaller time windows for large files
- Close other applications to free up memory

## API Costs

Estimated costs for AI features (using GPT-3.5-turbo):

- Per query: ~$0.0005 (0.05 cents)
- 100 queries: ~$0.05
- 1000 queries: ~$0.50

OpenAI provides $5-18 in free credits for new accounts.

## Development

### Adding New Features

1. Create new UI file in `ui/`
2. Create corresponding server file in `server/`
3. Source files in `app.R`
4. Update navigation in `ui_main.R`

### Testing

Test the AI integration:

```r
source("shiny_app/utils/ai_integration.R")
response <- get_ai_suggestion("What is ASAP?")
print(response)
```

## Support

- **ASAP Package**: https://github.com/LXiao06/ASAP
- **Issues**: Report bugs on GitHub
- **Documentation**: See ASAP package vignettes

## License

This Shiny app is part of the ASAP package and follows the same license.

## Acknowledgments

Built with:
- [Shiny](https://shiny.rstudio.com/) - Web framework
- [bslib](https://rstudio.github.io/bslib/) - Bootstrap 5 theming
- [OpenAI API](https://openai.com/) - AI assistance
- [ASAP](https://github.com/LXiao06/ASAP) - Audio analysis package
