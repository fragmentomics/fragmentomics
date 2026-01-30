# FragMentor Dockerfile
# Multi-stage build for minimal production image

# ============== Build Stage ==============
FROM python:3.12-slim AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Create virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python dependencies
COPY pyproject.toml .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir build && \
    pip install --no-cache-dir pysam polars numpy scipy matplotlib click rich

# Copy source and install
COPY src/ src/
COPY README.md LICENSE ./
RUN pip install --no-cache-dir .

# ============== Production Stage ==============
FROM python:3.12-slim AS production

LABEL maintainer="FragMentor Team"
LABEL description="The definitive toolkit for cfDNA fragmentomics analysis"
LABEL version="0.1.0"

# Install runtime dependencies only
RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Create non-root user
RUN useradd --create-home --shell /bin/bash fraguser
USER fraguser
WORKDIR /home/fraguser

# Default command
ENTRYPOINT ["fragmentomics"]
CMD ["--help"]

# ============== Development Stage ==============
FROM builder AS development

# Install dev dependencies
RUN pip install --no-cache-dir pytest pytest-cov ruff mypy

# Copy full source including tests
COPY . /app
WORKDIR /app

# Run tests by default in dev mode
CMD ["pytest", "-v"]
