FROM ghcr.io/mathworks-ref-arch/matlab-dockerfile/matlab:r2024a
USER root

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/
ENV UV_COMPILE_BYTECODE=1 UV_LINK_MODE=copy
ENV UV_PROJECT_ENVIRONMENT=/rte_dataset_venv

WORKDIR /app/rte_dataset
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    --mount=type=bind,source=.python-version,target=.python-version \
    uv sync --no-install-project --no-dev --frozen

ENV PATH="/rte_dataset_venv/bin:$PATH"
