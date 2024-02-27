# hash:sha256:3499650d8cdc294a7ef6c6bf6676c63c4836fd16413ca7d527531da5131fb76c
FROM registry.codeocean.allenneuraldynamics.org/codeocean/c2-seurat:latest

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
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
