#!/bin/bash
#PBS -N run_model
#PBS -l select=1:ncpus=1:mem=8gb:scratch_local=50gb
#PBS -l walltime=1:00:00 

#clean scratch
trap "clean_scratch" EXIT

#setting up variables
DATADIR=/auto/brno2/home/lojzero/MatabolismModel_Nlan_peroxisome_update
CONDA_ENV=/storage/brno2/home/lojzero/condaEnvs
CONDA_ENV2=/storage/plzen1/home/lojzero/condaEnvs
export TMPDIR=$SCRATCHDIR #changing scratchdir tmp
export GRB_LICENSE_FILE=$CONDA_ENV/iNlan20/gurobi_licence/gurobi.lic #you need to get your own license and put a path to it here!!!!!!!!

#go to the SCRATCHDIR
cd $SCRATCHDIR

#clone the iNlan20P repository
git clone https://github.com/LojzaZ/iNlan20P.git

#cd into the script dir
cd $SCRATCHDIR/iNlan20P/iNlan20P/ReconstructionAndAnalysis

#clone the original model from the iNlan20 repository
wget "https://raw.githubusercontent.com/stelmo/iNlan20/c13f06c073db560ccb46062ed26b46eb80ff63de/ReconstructionAndAnalysis/iNlan20.xml"

#activate env
module add mambaforge
mamba activate $CONDA_ENV/iNlan20

#make the results again
jupyter nbconvert \
    --to notebook \
    --execute \
    "Reconstruction of N. lanati GEM.ipynb"

#make the model
jupyter nbconvert \
    --to notebook \
    --execute \
    "Write model.ipynb"

#run the models
python final_script_run_models.py

#copy the results of the original model into its directory
mkdir -p res_iNlan20
cp *original* res_iNlan20

#copy the results of the updated model into its own directory
mkdir -p res_iNlan20P
cp *updated* res_iNlan20P

#copy the results of the reduction of FAOX flux into its own directory
mkdir -p res_FAOX_stepwise_reduction
cp *FAOXx* res_FAOX_stepwise_reduction

#finish the job
mamba deactivate
clean_scratch
