# R Project: Performance Analysis 
# (data.table vs. dplyr vs. SQL)

Author: **Marta Paira**

This repository contains the complete analysis for the bioinformatics project, aiming to compare the performance of `data.table`, `dplyr` (`tidyverse`), and `sqldf` (SQL), suggested as alternative package. The project covers 12 data wrangling tasks and a final single-cell analysis revision.

The project is structured to include:

* An **R benchmark script** (`progetto_completo.R`).

* An **R Package** (`ProgettoBioinfo`) containing the optimized (`data.table`) functions.

* Roxygen2 documentation (in the `man/` folder).

* A reproducible **Docker** environment (`Dockerfile`).

* An **R Markdown Report** (`report_progetto.Rmd`) which analyzes and visualizes the performance results.

* The CSV output files for every task (e.g., `Task1_FilterSummarize_Results.csv`).


## 1. Project Structure

The root directory (`project_oct25/`) contains the following key elements:

* `Dockerfile`: Defines the Ubuntu + R + R dependencies environment.

* `progetto_completo.R`: The R script that runs the entire benchmark analysis (batch mode).

* `report_progetto.Rmd`: The final R Markdown report that analyzes and visualizes the results.

* `DESCRIPTION`, `NAMESPACE`: Metadata files for the R package.

* `R/`: Folder containing the pure package functions (e.g., `task_functions.R`).

* `man/`: Generated Roxygen2 documentation (the `.Rd` files).

* `*.csv` (Data Files): The raw 14 data files (e.g., `bulk_counts_long.csv`, etc.).

* `results/`: Folder containing all generated output files (`Task*.csv`and `FinalTask*.png`).

---

## 2. Requirements

To replicate this analysis, the only external dependency required on your local machine is **Docker Desktop**.

* **[Docker Desktop](https://www.docker.com/products/docker-desktop/)**: This project is fully containerized. All R dependencies are installed and managed by the `Dockerfile`. 

* **Data Files:** The analysis requires the 14 `.csv` data files to be present in the root directory.

Docker is used to ensure a consistent and reproducible environment for the analysis. It allows us to create a virtual environment that contains all the necessary software and dependencies (R, `data.table`, `ggplot2`, etc.). 
This approach avoids issues related to different package versions or operating systems, ensuring the analysis runs smoothly and produces consistent results, regardless of the user's local machine configuration.
The basis for this environment is the **Dockerfile**, which serves as the blueprint for building the Docker image.

---

## 3. How to Run the Analysis 

This project is designed to be run entirely within the Docker container. You must have Docker installed and running.

### Step 1: Build the Docker Image

First, you must build the Docker image. This command reads the `Dockerfile`, downloads the Ubuntu base image, installs R, and installs all required R packages. This step only needs to be done once.

Open your terminal in the project's root directory and run:

```bash
docker build -t bioinfo-progetto-r:latest .
```

### Step 2: Run the Full Analysis Script (Generates Data)
This is the main command. It runs the entire `progetto_completo.R` script from start to finish. It will generate all result CSVs and the final plot (.png) and save them to the **`results/`** folder.

This command mounts your current project directory (represented by `$(pwd)`) into the `/app` directory *inside* the container. This is how the script finds your data and how the results are saved back to your computer.

```bash
# For macOS/Linux:
docker run --rm -v "$(pwd)":/app bioinfo-progetto-r:latest Rscript /app/progetto_completo.R
    
# For Windows (Command Prompt):
docker run --rm -v "%cd%":/app bioinfo-progetto-r:latest Rscript /app/progetto_completo.R
```

* `--rm`: This is a cleanup command. It automatically deletes the container when the script finishes.

* `-v "$(pwd)":/app`: This syncs your local folder with the container's working directory.

* You will see the `cat()` outputs from the R script in your terminal as it runs each Task. When the script is finished, all result files (`Task*.csv`, `FinalTask_CellTypeDistribution.png`) will be present in your local project folder.

---

## 4. R Package Details (ProgettoBioinfo)

The core logic of this analysis is encapsulated within a custom R package named ProgettoBioinfo.

You do not need to install or build this package manually. The `Dockerfile` handles the installation of this package inside the container automatically when the image is built.

The package serves two purposes:

* **Documentation**: It uses the Roxygen2 framework (files in `man/`) to document every optimized function created for this project.

* **Code Availability**: It ensures that the functions used for the `data.table` benchmarks are packaged as reusable components.

The package's source code is located in the `R/` folder, and the Roxygen2-generated documentation is in the `man/` folder.

Key functions provided by this package include:

* `task5_classify_lab_status()`: Performs the Non-Equi Join (interval join) for lab classification.

* `task6_rolling_join_vitals()`: Performs the nearest-time matching (Rolling Join).

* `task10_calculate_overlap_stats()`: Performs the foverlaps genomic spatial join.

* `finaltask_normalize_by_tissue()`: Calculates the normalized percentages for the single-cell analysis.
























