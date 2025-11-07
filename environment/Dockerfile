# hash:sha256:b02abc5aeaa30c8f518b73947791cd29fddedef3b13bb1cbde0889b37ce47606
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/c2-seurat:latest

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-askpass /

RUN mamba install -y \
        r-base==4.3.0 \
    && mamba clean -ya

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BiocVersion", \
        "IRanges", \
        "scater" \
    ))' # Original versions: 3.18.1 2.36.0 1.30.1

COPY postInstall /
RUN /postInstall
