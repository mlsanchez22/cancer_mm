FROM mambaorg/micromamba:1.5.3

# copy conda environment yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER ./environments/environment_py.yml /tmp/env.yaml

# allows using micromamba inside the docker
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# create conda environment
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# copy scripts 
#COPY --chown=$MAMBA_USER:$MAMBA_USER ./src/*.py /sbin

# Add conda PATH. Needed when singuarity execute the docker with "exec" instead of "run" command. (seen here: https://github.com/mamba-org/micromamba-docker/issues/155)
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
