FROM rocker/verse:4.2.1
LABEL maintainer="lorentzb"
RUN apt-get update -y && apt-get install -y  libcairo2-dev  libpng-dev  libjpeg-dev  zlib1g-dev  libicu-dev  pandoc  make  libcurl4-openssl-dev  libssl-dev  libxml2-dev  git  libfreetype6-dev  libfribidi-dev  libharfbuzz-dev  libfontconfig1-dev  libtiff-dev  libgit2-dev  libglpk-dev  libgmp3-dev && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(renv.config.pak.enabled = TRUE, repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages(c("renv","remotes","devtools"))'
RUN R -e 'devtools::install_github("kylebittinger/qiimer@7d7721315fbdfde808cecbe10661158ceec25e8b",build = TRUE)'
RUN R -e 'devtools::install_github("joey711/biom@d37425bee2c61189e6cfe24ff2c7e629966000de", build=TRUE)'
RUN R -e 'devtools::install_github("nick-youngblut/Tax4Fun@1.0",build=TRUE)'
RUN R -e 'install.packages("ppcor",INSTALL_opts="--no-lock")'
COPY renv.lock renv.lock
RUN R -e 'renv::restore()'
#CMD ["bash"]