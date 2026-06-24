#!/bin/bash
#PBS -N iNlan20_conda_env_setup
#PBS -l select=1:ncpus=1:mem=16gb:scratch_local=15gb
#PBS -l walltime=3:00:00

#trap clean scratch
trap "clean_scratch" EXIT

#path to your conda enviroments
CONDA_ENV=/storage/brno2/home/lojzero/condaEnvs

#set tmp and caches
export TMPDIR=$SCRATCHDIR
export CONDA_PKGS_DIRS=$SCRATCHDIR/conda_pkgs
export PIP_CACHE_DIR=$SCRATCHDIR/pip_cache
export XDG_CACHE_HOME=$SCRATCHDIR/xdg_cache
export MAMBA_ROOT_PREFIX=$SCRATCHDIR/mamba_root

#go into scratch
cd $SCRATCHDIR

#set env name and python and julia versions
ENV_NAME="$CONDA_ENV/iNlan20"
PYTHON_VERSION="3.8.3"
#JULIA_VERSION="1.5" #musi byt vic nez 1.5

#load mamba
module add mambaforge

#create env
mamba create -p "$ENV_NAME" -y -c anaconda python="$PYTHON_VERSION"

#install the conda env
mamba activate $ENV_NAME
mamba install -y -c conda-forge -c bioconda "xlrd=1.2.0" \
   cobra \
   pandas \
   numpy \
   openpyxl \
   "xlrd=1.2.0" \
   jupyter \
   jupyterlab \
   ipykernel \
   optlang


# Gurobi requires a valid licence — gurobipy is the Python interface = for free for stundents
pip install gurobipy

# NOTE: gurobipy is installed above but a Gurobi licence file (gurobi.lic)
# must be present on your system (usually at ~/gurobi.lic).
# Academic licences are free at https://www.gurobi.com/academia/academic-program-and-licenses/

#install julia and its packages packages
#mamba install -y -c conda-forge julia="$JULIA_VERSION"
mamba install -y -c conda-forge julia

julia --project=@. -e '
import Pkg
Pkg.add([
    "JSON",
    "ExcelReaders",
    "DataValues",
    "BioSequences",
    "DataFrames",
    "StatsBase",
    "FASTX"
])
Pkg.instantiate()
println("Julia packages installed successfully.")
'

# Register the Julia kernel so Jupyter can use it
julia -e '
import Pkg
Pkg.add("IJulia")
'

kernel_path=$(julia -e '
using IJulia
installkernel("julia-1.5")
' 2>&1 | grep -oP '(?<=kernelspec in ).*')


#check the kernel before moving
jupyter kernelspec list

#move the kernel to the conda env
cp -r "$kernel_path" $ENV_NAME/share/jupyter/kernels/julia-1.5

#check the kernel after moving
jupyter kernelspec list

##Registering Python kernel '$ENV_NAME' with Jupyter
ENV_NAME_SHORT=$(basename $ENV_NAME)

python -m ipykernel install \
    --user \
    --name "$ENV_NAME_SHORT" \
    --display-name "Python ($ENV_NAME_SHORT)"

ENV_NAME_SHORT="${ENV_NAME_SHORT,,}"

kernel_path=$(jupyter kernelspec list | grep inlan20 | awk '{print $2}')

cp -r "$kernel_path" $ENV_NAME/share/jupyter/kernels/inlan20

mamba deactivate
clean_scratch
