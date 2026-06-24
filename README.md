# iNlan20P
This repository contains the updated metabolic model of *Neocallimastix lanati*. Main update was including the compartment of peroxisomes into the model. See the paper with the original model iNlan20 (https://msystems.asm.org/content/6/1/e00002-21) and the paper with the updated model iNlan20P (add link)

# To run the model
1. Install the necessary conda environment (install_conda.sh)
2. Get your gurobi license at https://www.gurobi.com/academics
3. Use the script iNlan20P/ReconstructionAndAnalysis/run_model_iNlan20P to create the model and to get the data needed to produce the graphs in add link
4. Precomputed data are found in iNlan20P/ReconstructionAndAnalysis/res*
5. Precomputed model is located at iNlan20P/ReconstructionAndAnalysis/iNlan20P.xls and iNlan20P/ReconstructionAndAnalysis/iNlan20P.json

# Note
The original model was run on metacentrum (https://metavo.metacentrum.cz/) with PBS submission system. For your needs, you may need to change some paths (make sure the path to the gurobi licence is correct)
