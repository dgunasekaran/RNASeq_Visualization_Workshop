# RNA-seq QC & Visualization Workshop

Installation and setup for the workshop scripts and the RNA-seq QC Explorer app. Finish these steps **before** the workshop: they take 30–45 minutes, mostly waiting for packages to install.

You need a laptop with at least 8 GB RAM and 3 GB free disk space.

## 1. Install R (4.6.x)

Download R from CRAN and install it **before** RStudio:

| System | Download |
|---|---|
| Windows 10/11 | [CRAN › Download R for Windows › base](https://cloud.r-project.org/bin/windows/base/) |
| macOS, Apple silicon (M1–M4) | [CRAN › macOS](https://cloud.r-project.org/bin/macosx/), the **arm64** `.pkg` |
| macOS, Intel | [CRAN › macOS](https://cloud.r-project.org/bin/macosx/), the **x86_64** `.pkg` |
| Ubuntu / Debian | [CRAN Ubuntu instructions](https://cloud.r-project.org/bin/linux/ubuntu/) |

If you already have R 4.5 or older, update: Bioconductor packages are built for one R version at a time.

## 2. Install RStudio Desktop

1. Open the [RStudio download page](https://docs.posit.co/ide/user/download.html) and, under **RStudio Desktop**, pick the installer for your system (not RStudio Server or Pro).
2. **Windows:** run the `.exe` and accept the defaults.
   **macOS:** open the `.dmg`, drag **RStudio** into **Applications**, then eject the disk image.
   **Ubuntu / Debian:** `sudo apt install ./rstudio-*.deb` · **Fedora / RHEL:** `sudo dnf install ./rstudio-*.rpm`
3. Open RStudio (not the plain R app). The Console should start with `R version 4.6.x`.
4. Recommended: *Tools › Global Options › General*: untick "Restore .RData into workspace at startup" and set "Save workspace to .RData on exit" to **Never**.

Current RStudio supports Windows 11, macOS 14 or newer, and recent Ubuntu, Debian, Fedora and RHEL.

## 3. Get this repository

- **Without git:** click **Code › Download ZIP** at the top of this page and unzip it somewhere without spaces or accents in the path (for example `Documents/RNASeq_Visualization_Workshop`).
- **With git:** in RStudio choose *File › New Project › Version Control › Git* and paste this repository's URL.

Then open the project by double-clicking **`RNASeq_Visualization_Workshop.Rproj`**. Every script assumes this folder is the working directory.

## 4. Install the R packages

With the project open, run this once in the RStudio Console (10–25 minutes on a fresh R):

```r
source("scripts/install_workshop.R")
```

- If R asks *"Update all/some/none?"*, type `n`.
- If it asks *"Do you want to install from sources the packages which need compilation?"*, type `no`.

You are done when the Console prints **"Workshop packages installed …"**, then **"Loaded 7126 features x 24 samples … You're ready!"**, and a bar chart of library sizes appears in the Plots pane. A line listing missing *optional* packages is fine.

`install.R` installs only the app's packages. `scripts/install_workshop.R` runs it and adds the packages the workshop scripts need.

## 5. Check the data and the app

In the RStudio Console, with the project open:

```r
source("scripts/workshop_helpers.R")
counts  <- read_counts("data/demo_counts.csv")
coldata <- read_metadata("data/demo_metadata.csv", colnames(counts))
dim(counts)                   # 7126   24
table(coldata$condition)      # 4 samples in each of 6 groups
shiny::runApp()               # opens the app; first load takes about a minute
```

In the app:

1. The **Overview** tab shows 24 samples and 7,126 features.
2. Type **HSP104** in the sidebar gene search box. It should be highlighted on the volcano plot (**Differential expression** tab) and plotted in the **Gene explorer** tab.
3. Close the app with the red **Stop** button in the RStudio Console.

## Troubleshooting

| Symptom | Fix |
|---|---|
| `cannot open file 'R/helpers.R'` or "Run the workshop scripts from the repo folder" | Open `RNASeq_Visualization_Workshop.Rproj`, or *Session › Set Working Directory › To Project Directory* |
| `Bioconductor version '3.xx' requires R version '4.x'` | Install R 4.6.x, restart RStudio, re-run `source("scripts/install_workshop.R")` |
| `package 'X' is not available for this version of R` | Check `R.version.string`, then `BiocManager::install("X")` |
| `lazy-load database ... is corrupt` | *Session › Restart R*, then run again |
| Windows: `cannot remove prior installation` / permission denied | Close other R sessions; accept the personal library R offers |
| macOS: `xcrun: error` or compiler errors | Answer `no` to compiling from source; check you installed the right arm64 / x86_64 R |
| Linux: compilation of `xml2`, `curl` or `ragg` fails | `sudo apt install libxml2-dev libcurl4-openssl-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libpng-dev` |
| App opens but plots are blank | Click *Open in Browser* at the top of the app window |

Still stuck? Send the output of `sessionInfo()` and the full error message to the organizer, Deepika Gunasekaran (dgunasekaran@ucmerced.edu).
