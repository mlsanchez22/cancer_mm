FROM mambaorg/micromamba:1.5.3

# copy conda environment
COPY --chown=$MAMBA_USER:$MAMBA_USER ./environments/environment_r.yml /tmp/env.yaml

# create conda environment
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# copy scripts 
#COPY --chown=$MAMBA_USER:$MAMBA_USER ./src/*.R /sbin

# Add conda PATH. Needed when singuarity execute the docker with "exec" instead of "run" command. (seen here: https://github.com/mamba-org/micromamba-docker/issues/155)
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
