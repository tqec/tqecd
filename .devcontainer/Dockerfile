FROM mcr.microsoft.com/devcontainers/python:3.12

# Installing pandoc, required for documentation generation.
RUN <<EOF
apt-get update
apt-get install -y --no-install-recommends pandoc graphviz
EOF
