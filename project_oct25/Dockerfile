# 1. SET THE BASE IMAGE
# Defines the base image upon which subsequent layers (R and dependencies) will be installed. 
# We use ubuntu:latest (Noble Numbat 24.04), the most recent base.

FROM ubuntu:latest


# SET THE ENVIRONMENT VARIABLE FOR NON-INTERACTIVITY
# To prevent commands (e.g., apt-get) from blocking while waiting for user input.

ENV DEBIAN_FRONTEND=noninteractive


# 2. CONFIGURE THE R CRAN REPOSITORY AND SECURITY KEYS (For the latest base)
# Adds the official R-CRAN repository to Ubuntu's sources, necessary to install the latest version of R (r-base).
# Updates the list of known packages in Ubuntu.
# Installs necessary tools to add new repositories (wget, gpg).
# Imports the R-CRAN security key (guarantee seal).
# Adds the R repository address (for Ubuntu Noble) to the apt list.

RUN apt-get update && \
    apt-get install -y --no-install-recommends dirmngr wget software-properties-common && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor | tee /usr/share/keyrings/cran.gpg > /dev/null && \
    # Add CRAN repository for Noble (latest version)
    echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" | tee /etc/apt/sources.list.d/cran.list


# 3. INSTALL R, COMPILERS, AND ESSENTIAL DEPENDENCIES
# Updates the package list.
# Installs the software in non-interactive (-y) and minimal mode (do not install recommended packages).
# pandoc: essential for R Markdown (to convert .Rmd to .html/pdf).
# Installs R, build-essential (e.g., C/C++ compilers), dependencies, and required packages.
# Cleans the 'apt' cache to reduce the final image size.

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        curl grep gdebi-core pandoc \
        r-base r-base-dev build-essential libcurl4-gnutls-dev libssl-dev libxml2-dev libpq-dev sudo psmisc \
        libsqlite3-dev \ 
        && rm -rf /var/lib/apt/lists/*


# 4. INSTALL R PACKAGES FOR THE PROJECT
# Sets the default working directory inside the container.
# Executes R and installs the necessary R packages for the analysis.

WORKDIR /app
RUN R -e "install.packages(c('data.table', 'dplyr', 'tidyverse', 'knitr', 'rmarkdown', 'sqldf', 'tidyr', 'ggplot2', 'scales'), repos = 'https://cloud.r-project.org/')"


# 4.5. INSTALL CUSTOM R PACKAGE (ProgettoBioinfo)
# The source package files must be present inside the image to be installed.
# This installs the package so the final script can use library(ProgettoBioinfo).
# We copy the essential files (R/, DESCRIPTION, NAMESPACE) from the local machine into the Docker image's working directory (/app).
COPY R /app/R
COPY DESCRIPTION /app/DESCRIPTION
COPY NAMESPACE /app/NAMESPACE

# Use R CMD INSTALL for reliable package installation from local source.
RUN R CMD INSTALL /app


# 5. CONFIGURE USER AND WORKDIR
# Creates a non-root 'rstudio' user for security reasons.
# Sets the password for the user ('user:password' format).
# Sets the user's default shell to /bin/bash (more convenient than /bin/sh).

RUN useradd -m rstudio && \
    echo "rstudio:rstudio" | chpasswd && \
    usermod -s /bin/bash rstudio


# 6. EXPOSE THE PORT and set the default command to BASH
# Informs Docker that the container uses port 8787 (RStudio standard).
# Sets /bin/bash (terminal) as the default command when starting the container.

EXPOSE 8787
CMD ["/bin/bash"]


## docker build -t bioinfo-progetto-r:latest .

## docker run --rm -v /Users/martapaira/Desktop/project_oct25:/app bioinfo-progetto-r:latest Rscript /app/progetto_completo.R