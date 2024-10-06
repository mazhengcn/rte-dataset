FROM ghcr.io/mathworks-ref-arch/matlab-dockerfile/matlab:r2024a

COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

ENV PYTHONUNBUFFERED=True UV_LINK_MODE=copy
ENV UV_PROJECT_ENVIRONMENT=/opt/rte-dataset

USER root

RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    --mount=type=bind,source=.python-version,target=.python-version \
    uv sync --no-install-project --no-dev

ENV PATH="/opt/rte-dataset/bin:$PATH"
