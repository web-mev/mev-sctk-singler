FROM ghcr.io/web-mev/base-sctk-docker:sha-b6f556d679aaaa305a80a725e3813431fb93d12f

# Install the optparse package, which allows better commandline arg parsing
RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN mkdir -p /opt/software
ADD singler.R /opt/software/
ADD create_cache.R /opt/software/

# Caches the databases so we don't depend on live downloads during runtime.
Run Rscript /opt/software/create_cache.R