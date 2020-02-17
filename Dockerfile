## Start from this Docker image
FROM r-base:3.6.2
#FROM ubuntu

## Install R in Docker image
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y locales \
    && sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && dpkg-reconfigure --frontend=noninteractive locales \
    && update-locale LANG=en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'
RUN echo en_US.UTF-8 UTF-8 >> /etc/locale.gen && locale-gen
ENV PYTHONIOENCODING=utf-8

#RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y r-base

## Install R packages in Docker image
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install("S4Vectors")'
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('readr')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('reshape2')"
RUN Rscript -e "install.packages('caret')"
RUN Rscript -e "install.packages('e1071')"
RUN Rscript -e "install.packages('randomForest')"
RUN Rscript -e "install.packages('import')"
RUN Rscript -e "install.packages('foreach')"
RUN Rscript -e "install.packages('xgboost')"

## Copy your files into Docker image
USER root
RUN mkdir /rds
COPY model-fits.rds /rds/
COPY dnaseqGenes.rds /rds/
COPY rnaseqGenes.rds /rds/
COPY rna-preProc.rds /rds/
COPY beatAML_pred.R /usr/local/bin/
RUN chmod a+x /usr/local/bin/beatAML_pred.R

## Make Docker container executable
ENTRYPOINT Rscript --vanilla /usr/local/bin/beatAML_pred.R
