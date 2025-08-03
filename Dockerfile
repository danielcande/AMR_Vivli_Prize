# Use the official R base image
FROM r-base:4.3.0

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgdal-dev \
    libudunits2-dev \
    libproj-dev \
    libgeos-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libxt-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libfreetype-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    pandoc \
    pkg-config \
    cmake \
    make \
    g++ \
    gcc \
    libc6-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages step by step to handle dependencies better
RUN R -e "install.packages('remotes', repos='https://cran.rstudio.com/')"

# Install core packages first
RUN R -e "install.packages(c('shiny', 'shinythemes'), repos='https://cran.rstudio.com/')"

# Install system font packages with specific configuration
RUN R -e "install.packages('systemfonts', repos='https://cran.rstudio.com/', configure.args='--with-freetype-includes=/usr/include/freetype2')"

# Install tidyverse and related packages
RUN R -e "install.packages(c('dplyr', 'ggplot2', 'tidyr', 'readr', 'purrr', 'tibble', 'stringr', 'forcats'), repos='https://cran.rstudio.com/')"

# Install spatial packages
RUN R -e "install.packages(c('sf', 'rnaturalearth'), repos='https://cran.rstudio.com/')"

# Install remaining packages
RUN R -e "install.packages(c('countrycode', 'xgboost', 'gt', 'nnet', 'scales', 'viridis', 'RColorBrewer'), repos='https://cran.rstudio.com/')"

# Copy the Shiny app files
COPY app.R /app/
COPY Vivli\(corrected\).RData /app/

# Create a non-root user
RUN useradd -m -u 1000 user
RUN chown -R user:user /app
USER user

# Expose the port that the app will run on
EXPOSE 7860

# Set environment variables for Shiny
ENV SHINY_HOST=0.0.0.0
ENV SHINY_PORT=7860

# Start the Shiny app directly
CMD ["R", "-e", "options(shiny.host='0.0.0.0', shiny.port=7860); source('app.R')"]