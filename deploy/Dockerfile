# CIFcurator server image. Built by docker-compose.yml from the repository root:
#   docker compose -f deploy/docker-compose.yml up -d --build
# Python packages, MMseqs2 and MUSCLE come from conda (environment-server.yml); the CIFcurator
# code (scripts/, src/, server/) is copied in at build time, so rebuild after changing it.
FROM mambaorg/micromamba:2.10.0-debian12-slim

COPY --chown=$MAMBA_USER:$MAMBA_USER deploy/environment-server.yml /tmp/environment-server.yml
RUN micromamba install -y -n base -f /tmp/environment-server.yml && micromamba clean --all --yes

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/ scripts/
COPY --chown=$MAMBA_USER:$MAMBA_USER src/ src/
COPY --chown=$MAMBA_USER:$MAMBA_USER server/ server/

# Jobs and the PDB sequence file live in /data (a Docker volume, kept when the image is rebuilt)
USER root
RUN mkdir -p /data && chown $MAMBA_USER:$MAMBA_USER /data && rm -rf server/server_data server/__pycache__
USER $MAMBA_USER

# Commit of the code, written to every job's logs/parameters.json (there is no .git in the image)
ARG CIFCURATOR_COMMIT=unknown
ENV CIFCURATOR_COMMIT=$CIFCURATOR_COMMIT PYTHONUNBUFFERED=1

EXPOSE 8000
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "server/cifcurator_server.py", "--host", "0.0.0.0", "--port", "8000", "--data-dir", "/data"]
