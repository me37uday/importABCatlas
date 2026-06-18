FROM rocker/tidyverse

# install reticulate
RUN R -e "install.packages('reticulate', repos='https://cloud.r-project.org/')"

# install seurat
RUN R -e "install.packages('Seurat', repos='https://cloud.r-project.org/')"

# install python 3.10
RUN apt-get update && apt-get install -y python3.10 python3-venv python3-pip python3-dev


COPY . /app

WORKDIR /app

#RUN R -e "devtools::load_all('.')"

#CMD ["R", "-e", "whb <- load_data_WHB()"]