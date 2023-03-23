FROM rocker/verse:4.2.0
LABEL maintainer="lorentzb"
RUN apt-get update -y && apt-get install -y libicu-dev pandoc make zlib1g-dev libcurl4-openssl-dev  git  libssl-dev  libfreetype6-dev  libfribidi-dev  libharfbuzz-dev  libxml2-dev  libfontconfig1-dev  libgit2-dev  libjpeg-dev  libpng-dev libtiff-dev  libglpk-dev  libgmp3-dev && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(renv.config.pak.enabled = TRUE, repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages(c("renv","remotes"))'
COPY renv.lock renv.lock
RUN R -e 'renv::restore()'
CMD [ 'bash' ]
