#!/usr/bin/env bash
set -euo pipefail

# Bootstrap the driver environment for Ubuntu under WSL2.
# The individual Snakemake rules create their own pinned Conda environments.

if [[ "$(uname -s)" != "Linux" ]]; then
  echo "This setup script is intended for Linux/WSL2." >&2
  exit 1
fi

ARCH="$(uname -m)"
if [[ "$ARCH" != "x86_64" ]]; then
  echo "WARNING: Minigraph-Cactus is currently intended for x86-64. Detected: $ARCH" >&2
fi

sudo apt-get update
sudo apt-get install -y \
  git curl wget ca-certificates bzip2 unzip build-essential

MINIFORGE="$HOME/miniforge3"
if ! command -v conda >/dev/null 2>&1; then
  if [[ ! -x "$MINIFORGE/bin/conda" ]]; then
    echo "Installing Miniforge to $MINIFORGE"
    tmp="$(mktemp -d)"
    curl -L -o "$tmp/miniforge.sh" \
      https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash "$tmp/miniforge.sh" -b -p "$MINIFORGE"
    rm -rf "$tmp"
  fi
  export PATH="$MINIFORGE/bin:$PATH"
fi

# Make conda available in this non-interactive shell.
CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda config --set channel_priority strict
conda create -y -n pangenome-hubo \
  -c conda-forge -c bioconda \
  python=3.11 snakemake

conda activate pangenome-hubo

echo
echo "Setup complete."
echo "Snakemake: $(snakemake --version)"
echo
echo "Next commands:"
echo "  conda activate pangenome-hubo"
echo "  snakemake -n -p --software-deployment-method conda --cores 4"
echo "  snakemake -p --software-deployment-method conda --cores 8"
echo
echo "Before the real run, put FASTA assemblies in data/assemblies/ and edit samples.tsv/config.yaml."
