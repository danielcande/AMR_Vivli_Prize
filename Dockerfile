# Use rocker/r-ver for a cleaner base
FROM rocker/r-ver:4.3.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgdal-dev \
    libudunits2-dev \
    libproj-dev \
    libgeos-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages in stages
RUN R -e "install.packages('shiny', repos='https://cran.rstudio.com/')"
RUN R -e "install.packages('tidyverse', repos='https://cran.rstudio.com/')"
RUN R -e "install.packages('shinythemes', repos='https://cran.rstudio.com/')"
RUN R -e "install.packages('sf', repos='https://cran.rstudio.com/')"
RUN R -e "install.packages('countrycode', repos='https://cran.rstudio.com/')"
RUN R -e "install.packages('gt', repos='https://cran.rstudio.com/')"
RUN R -e "install.packages('nnet', repos='https://cran.rstudio.com/')"
RUN R -e "install.packages('scales', repos='https://cran.rstudio.com/')"

# Install rnaturalearth without problematic dependencies
RUN R -e "install.packages('rnaturalearth', repos='https://cran.rstudio.com/', dependencies=c('Depends', 'Imports'))"

# Install XGBoost
RUN R -e "install.packages('xgboost', repos='https://cran.rstudio.com/')"

# Set working directory
WORKDIR /app

# Copy your app files
COPY . /app/

# Expose port
EXPOSE 3838

# Run the Shiny app directly
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]