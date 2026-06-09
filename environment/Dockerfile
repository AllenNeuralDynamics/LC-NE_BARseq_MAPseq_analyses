# hash:sha256:ed427fb9cb4bf670ddcf31adbaf6cc4a91eec73fdd14d9fe8dd3ed6ab7ef62ea
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/c2-seurat:latest

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

COPY barseq-r4.yml /tmp/barseq-r4.yml

RUN mamba update -n base -c conda-forge mamba libsolv -y || true \
 && rm -rf /opt/conda/envs/* /etc/mamba-envs/* \
 && mamba env create -f /tmp/barseq-r4.yml \
 && rm /tmp/barseq-r4.yml \
 && echo "/opt/conda/envs/barseq-r4/lib" > /etc/ld.so.conf.d/barseq-r4.conf \
 && ldconfig

# aind-data-schema (Python) is used by code/08_generate_metadata.py to write
# AIND-compliant derived-asset metadata into /results/.
RUN python3 -m pip install --no-cache-dir aind-data-schema==2.7.1 aind-data-schema-models==5.4.1
