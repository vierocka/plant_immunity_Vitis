#!/bin/bash -e
# AIM: install a local, standalone ColabFold (AlphaFold2-Multimer) environment.
# Based on YoshitakaMo/localcolabfold's public installer
# (https://github.com/YoshitakaMo/localcolabfold), used as-is.
# Run once, on a machine/HPC node with a GPU.

type wget 2>/dev/null || { echo "wget is not installed. Please install it using apt or yum." ; exit 1 ; }

CURRENTPATH=$(pwd)
COLABFOLDDIR="${CURRENTPATH}/localcolabfold"
mkdir -p "${COLABFOLDDIR}"
cd "${COLABFOLDDIR}"

wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash ./Miniforge3-Linux-x86_64.sh -b -p "${COLABFOLDDIR}/conda"
rm Miniforge3-Linux-x86_64.sh

source "${COLABFOLDDIR}/conda/etc/profile.d/conda.sh"
export PATH="${COLABFOLDDIR}/conda/condabin:${PATH}"
conda update -n base conda -y
conda create -p "$COLABFOLDDIR/colabfold-conda" -c conda-forge -c bioconda \
    git python=3.10 openmm==8.2.0 pdbfixer \
    kalign2=2.04 hhsuite=3.3.0 mmseqs2 -y
conda activate "$COLABFOLDDIR/colabfold-conda"

"$COLABFOLDDIR/colabfold-conda/bin/pip" install --no-warn-conflicts \
    "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"
"$COLABFOLDDIR/colabfold-conda/bin/pip" install "colabfold[alphafold]"
"$COLABFOLDDIR/colabfold-conda/bin/pip" install --upgrade "jax[cuda12]==0.5.3"

echo "Done. Activate with: conda activate ${COLABFOLDDIR}/colabfold-conda"
echo "Binaries: colabfold_batch, colabfold_search"
