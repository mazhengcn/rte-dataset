FROM ghcr.io/mathworks-ref-arch/matlab-dockerfile/matlab:r2024a

USER root
# Install system packages
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y --no-install-recommends git screen \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:latest --chown=matlab: /uv /uvx /bin/
ENV UV_LINK_MODE=copy
# Only use the managed Python version
ENV UV_PYTHON_PREFERENCE=only-managed
# Configure the Python directory so it is consistent
ENV UV_PYTHON_INSTALL_DIR=/python

WORKDIR /python
RUN --mount=type=bind,source=.python-version,target=.python-version \
    uv python install
