FROM mambaorg/micromamba:1.5.3

RUN micromamba install -y -n base -c conda-forge python=3.10 pip procps-ng && \
    micromamba clean --all --yes

RUN micromamba run -n base pip install drexml==1.1.1
# Add conda PATH. Needed when singuarity execute the docker with "exec" instead of "run" command. (seen here: https://github.com/mamba-org/micromamba-docker/issues/155)
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
