# Use Rocker Shiny image optimized for Shiny applications
FROM --platform=linux/amd64 rocker/shiny:4.4.0

# Disable renv to avoid conflicts with remotes
ENV RENV_ACTIVATE=FALSE
ENV RENV_PATHS_CACHE=OFF

# Install system dependencies in one layer with clean up
RUN apt-get update && apt-get install -y \
  libharfbuzz-dev \
  libfribidi-dev \
  libfontconfig1-dev \
  libfreetype6-dev \
  libpng-dev \
  libcairo2-dev \
  libtiff5-dev \
  libjpeg-dev \
  libxml2-dev \
  libssl-dev \
  libcurl4-openssl-dev \
  pandoc \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Expose port 3838 for Shiny
EXPOSE 3838

# Copy DESCRIPTION first for better caching
COPY DESCRIPTION ./

# Configure CRAN repository and install pak (faster than remotes)
RUN R --vanilla -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest')); \
    install.packages('pak', dependencies = TRUE)"

# Install project dependencies using pak (much faster)
RUN R --vanilla -e "pak::local_install_deps('.', dependencies = NA)"

# Copy only necessary application files
COPY R/ ./R/
COPY NAMESPACE ./

# Install the package
RUN R --vanilla -e "pak::local_install('.', dependencies=FALSE)"


COPY diann_dashboard/ ./diann_dashboard/

# Command to run the Shiny application
CMD ["R", "--vanilla", "-e", "library(diann); shiny::runApp('diann_dashboard/', port = 3838, host = '0.0.0.0')"]
