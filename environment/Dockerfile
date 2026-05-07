FROM registry.codeocean.allenneuraldynamics.org/codeocean/c2-seurat:latest
ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-ask-pass /
COPY barseq-r4.yml /etc/mamba-envs/barseq-r4.yml

# Build the barseq-r4 R conda environment.
# The `|| true` on the mamba update line lets the build continue if
# mamba's self-update fails -- the rest of the build doesn't depend on it.
RUN mamba update -n base -c conda-forge mamba libsolv -y || true \
 && rm -rf /opt/conda/envs/barseq-r4 \
 && mamba env create -f /etc/mamba-envs/barseq-r4.yml \
 && echo "/opt/conda/envs/barseq-r4/lib" > /etc/ld.so.conf.d/barseq-r4.conf \
 && ldconfig
