# TODO: try 1.5.6-focal-cuda-11.8.0 -> Creo que detecta los driver compatibles en el momento de compilar, y por eso no funciona en github-hosted actions
#FROM nvidia/cuda:11.8.0-devel-ubuntu20.04
FROM mambaorg/micromamba:1.5.3

RUN micromamba install -y -n base --override-channels \
  -c conda-forge \
  cuda cuda-nvcc cuda-toolkit gxx=11.2 python=3.10 pip procps-ng && \
  micromamba clean --all --yes

RUN micromamba run -n base pip install --no-cache-dir --no-binary=shap drexml==1.1.1

# Add conda PATH. Needed when singuarity execute the docker with "exec" instead of "run" command. (seen here: https://github.com/mamba-org/micromamba-docker/issues/155)
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"