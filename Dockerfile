FROM rocker/r-ver:4.3.3

# ── System libraries ──────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    git curl wget build-essential \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libgit2-dev libssh2-1-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    libgsl-dev libhdf5-dev \
    libbamtools-dev samtools \
    python3 python3-pip python3-venv \
    pandoc \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

# ── Python ────────────────────────────────────────────────────────────────────
RUN pip3 install --no-cache-dir jupyterlab numpy pandas && \
    pip3 install --no-cache-dir torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121 && \
    pip3 install --no-cache-dir pyro-ppl

# ── JupyterLab: no token, no password ────────────────────────────────────────
RUN jupyter lab --generate-config && \
    echo "c.ServerApp.open_browser = False"         >> /root/.jupyter/jupyter_lab_config.py && \
    echo "c.ServerApp.ip = '0.0.0.0'"              >> /root/.jupyter/jupyter_lab_config.py && \
    echo "c.ServerApp.root_dir = '/sharedFolder'"   >> /root/.jupyter/jupyter_lab_config.py && \
    echo "c.ServerApp.allow_root = True"            >> /root/.jupyter/jupyter_lab_config.py

# ── R: BiocManager + remotes (con stop() per rilevare errori) ─────────────────
RUN R -e "install.packages(c('BiocManager','remotes'), repos='https://cloud.r-project.org', Ncpus=parallel::detectCores()); \
          if (!requireNamespace('remotes', quietly=TRUE)) stop('remotes install failed')"

# IRkernel
RUN R -e "install.packages('IRkernel', repos='https://cloud.r-project.org'); \
          IRkernel::installspec(user=FALSE)"

# ── Bioconductor ──────────────────────────────────────────────────────────────
RUN R -e "BiocManager::install(c( \
    'GenomicAlignments','SummarizedExperiment','plyranges', \
    'Rsamtools','GenomeInfoDb', \
    'BSgenome.Hsapiens.UCSC.hg38', \
    'GenomicRanges','Biostrings', \
    'BiocGenerics','S4Vectors','GenomicFeatures', \
    'edgeR','ComplexHeatmap' \
  ), ask=FALSE, update=FALSE, Ncpus=parallel::detectCores())"

# ── CRAN ──────────────────────────────────────────────────────────────────────
RUN R -e "install.packages(c('ggdendro','ggplot2','dplyr','tidyr','patchwork'), \
            repos='https://cloud.r-project.org', Ncpus=parallel::detectCores())"

# ── GitHub via remotes ────────────────────────────────────────────────────────
RUN R -e "remotes::install_github('bowang-lab/simATAC',        upgrade='never')"
RUN R -e "remotes::install_github('Aelita-Stone/AtaCNV',       upgrade='never')"
RUN R -e "remotes::install_github('colomemaria/epiAneufinder', upgrade='never')"
RUN R -e "remotes::install_github('caravagnalab/rcongas',      upgrade='never')"

# ── Mount point ───────────────────────────────────────────────────────────────
RUN mkdir -p /sharedFolder
WORKDIR /sharedFolder
EXPOSE 8888

CMD ["sh", "-c", "jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root --ServerApp.token='' --ServerApp.password=\"$JUPYTER_PASSWORD_HASH\""]
