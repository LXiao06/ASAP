# Set up environment

## Step 1 — Set up environment

### Install R

ASAP requires **[R](https://cran.r-project.org) version 4.2 or higher**.
Download and install R for your operating system from the official CRAN
mirror.

### Install RStudio (Recommended IDE)

While ASAP can be used from any R environment, I recommend
**[RStudio](https://posit.co/download/rstudio-desktop)** as your
integrated development environment (IDE). It provides a comfortable
interface for running scripts, inspecting data, and viewing plots.

------------------------------------------------------------------------

## Step 2 — Install the `remotes` Package

ASAP is distributed through GitHub for now. The easiest way to install
GitHub-hosted R packages is via the `remotes` package.

Open RStudio and run the following in the Console:

``` r
install.packages("remotes")
```

------------------------------------------------------------------------

## Step 3 — Install ASAP

### Recommended: Latest Release

I recommend installing the latest **stable release** of ASAP for
everyday use. Releases are tagged and tested, so this is the safest
option.

``` r
remotes::install_github("LXiao06/ASAP@*release")
```

### Development Version

If you want access to the newest features before they are officially
released, you can install the **development version** directly from the
`main` branch. Note that development versions may contain experimental
features.

``` r
remotes::install_github("LXiao06/ASAP")
```

### A Specific Previous Version

If you need to reproduce results from an earlier version of ASAP, you
can install any past release by specifying its version tag (e.g.,
`v0.3.3`):

``` r
remotes::install_github("LXiao06/ASAP@v0.3.3")
```

You can browse all available version tags on the GitHub [Releases
page](https://github.com/LXiao06/ASAP/releases).

------------------------------------------------------------------------

## Step 4 — Load ASAP

Once installed, load the package in your R session:

``` r
library(ASAP)
```

------------------------------------------------------------------------

## Troubleshooting

**If installation fails**, a common culprit is missing system-level
dependencies required by audio packages (e.g., `tuneR`, `seewave`). Try
the following:

- **macOS**: Install Xcode command-line tools via
  `xcode-select --install` in Terminal.
- **Linux**: Install `libsndfile` via your system package manager (e.g.,
  `sudo apt-get install libsndfile1-dev`).
- **Windows**: Ensure
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is installed
  and on your PATH.

If you encounter a `package 'xxx' is not available` error, try updating
your R and `remotes` package first:

``` r
update.packages(ask = FALSE)
install.packages("remotes")
```

For further help, please open an issue on the [GitHub Issues
tracker](https://github.com/LXiao06/ASAP/issues).
