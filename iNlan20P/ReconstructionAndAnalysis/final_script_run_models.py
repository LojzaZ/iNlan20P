##############################
## load libraries           ##
##############################
import cobra
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
import json
import numpy as np
import gurobipy
from itertools import product

##############################
## configuration            ##
##############################
MODEL_PATHS = {
    "original" : "iNlan20.xml",
    "updated"  : "iNlan20P.xml",
    "updated_noFAOX"  : "iNlan20P.xml",
}

EXCHANGE_CONSTRAINTS = {
    "EX_h2_e"    : (0.047483008, 0.189248895),
    "EX_lac__D_e": (0.716237324, 1.08919967),
    "EX_for_e"   : (1.087540278, 1.793056312),
    "EX_etoh_e"  : (0.472808257, 1.013514567),
    "EX_ac_e"    : (0.423987759, 0.711876893),
    "EX_succ_e"  : (0.02143,     0.04529),
}

N_SAMPLES = 2000

##############################
## helper functions         ##
##############################
def load_flux_tibble(fbaloc, active_only=True, sort_by_magnitude=True):
    """Load FBA flux JSON into a tidy DataFrame."""
    with open(fbaloc, "r") as f:
        fluxes = json.load(f)
    flux_df = (
        pd.DataFrame.from_dict(fluxes, orient="index", columns=["flux"])
        .rename_axis("reaction_id")
        .reset_index()
    )
    if active_only:
        flux_df = flux_df[flux_df["flux"] != 0]
    if sort_by_magnitude:
        flux_df = flux_df.sort_values("flux", key=abs, ascending=False)
    return flux_df.reset_index(drop=True)

def get_model_reactions(model):
    """Extract all reactions and their reactants/products from a COBRApy model."""
    rows = []
    for rxn in model.reactions:
        for met, stoich in rxn.metabolites.items():
            rows.append({
                "reaction_id"  : rxn.id,
                "reaction_name": rxn.name,
                "metabolite_id": met.id,
                "metabolite_name": met.name,
                "stoichiometry": stoich,
                "role"         : "reactant" if stoich < 0 else "product"
            })
    return pd.DataFrame(rows)

def get_model_bounds(model):
    """Extract all reaction bounds from a COBRApy model."""
    return pd.DataFrame([
        {
            "reaction_id"  : rxn.id,
            "reaction_name": rxn.name,
            "lower_bound"  : rxn.lower_bound,
            "upper_bound"  : rxn.upper_bound,
        }
        for rxn in model.reactions
    ])

def apply_exchange_constraints(model, constraints, biomass = True):
    """Apply a dict of {rxn_id: (lb, ub)} exchange constraints to the model."""
    for rxn_id, (lb, ub) in constraints.items():
        rxn             = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = lb
        rxn.upper_bound = ub
    print(f"  Applied {len(constraints)} exchange constraints")

def apply_biomass_constrain(model, fba_sol, wiggle_room = 0.9):
    """Apply a biomass constrain so only the biologically relvant parts are explored during sampling."""
    mu = fba_sol.fluxes["Biomass"]
    bm_rxn = model.reactions.get_by_id('Biomass')
    mu_cons = model.problem.Constraint(bm_rxn.flux_expression, lb = wiggle_room*mu, ub=mu)
    model.add_cons_vars(mu_cons)
    print(f"  Applied biomass constrain")

def run_fba(model, tag):
    """Run FBA, save fluxes as JSON and TSV, print and save summary."""
    sol = model.optimize()
    print(f"  [{tag}] FBA status: {sol.status}  |  growth: {sol.fluxes['Biomass']:.6f}")
    if sol.status != "optimal":
        raise RuntimeError(f"FBA infeasible for {tag}")
    
    json_path = f"fluxes_{tag}.json"
    tsv_path  = f"fluxes_{tag}.tsv"
    sol.fluxes.to_json(json_path)
    load_flux_tibble(json_path).to_csv(tsv_path, sep="\t", index=False)
    
    summary = model.summary()
    print(summary)
    with open(f"summary_{tag}.txt", "w") as f:
        f.write(str(summary))
    
    return sol

def run_fva(model, tag, fraction=0.90):
    """Run FVA and save results as TSV."""
    fva = flux_variability_analysis(model, fraction_of_optimum=fraction)
    fva.to_csv(f"fva_{tag}.tsv", sep="\t")
    fixed    = (fva.maximum - fva.minimum).abs() < 1e-6
    print(f"  [{tag}] FVA done — {fixed.sum()} fixed / {(~fixed).sum()} variable reactions")
    return fva

def run_sampling(model, tag, n=N_SAMPLES):
    """Attempt ACHR sampling; warn and skip if warmup is empty."""
    sampler = cobra.sampling.ACHRSampler(model, thinning=100)
    print(f"  [{tag}] Warmup points: {sampler.n_warmup}")
    if sampler.n_warmup == 0:
        print(f"  [{tag}] WARNING — warmup empty, skipping sampling (constraints too tight)")
        return None
    s = cobra.sampling.sample(model, n, method="achr")
    s.to_json(f"samples_{tag}.json")
    print(f"  [{tag}] Sampling done — {n} samples saved")
    return s

def save_model_tables(model, prefix):
    """Save reaction/metabolite and bounds tables."""
    get_model_reactions(model).to_csv(f"{prefix}_reactions.tsv",  sep="\t", index=False)
    get_model_bounds(model).to_csv(   f"{prefix}_bounds.tsv",     sep="\t", index=False)
    print(f"  Saved reaction and bounds tables for {prefix}")

##############################
## model modifications      ##
##############################
def modify_updated_model(model):
    """All structural changes specific to the updated model."""
    
    # ── check which metabolites already exist ─────────────────────────────────
    existing = [m.id for m in model.metabolites]
    new_mets = {
        "hdec_e"  : Metabolite("hdec_e",   name="Hexadecanoic acid (extracellular)", compartment="e"),
        "hdec_c"  : Metabolite("hdec_c",   name="Hexadecanoic acid",                compartment="c"),
        "pmtcoa_c": Metabolite("pmtcoa_c", name="Palmitoyl-CoA",                    compartment="c"),
        "ppi_c"   : Metabolite("ppi_c",    name="Diphosphate",                      compartment="c"),
        "amp_c"   : Metabolite("amp_c",    name="AMP",                              compartment="c"),
    }
    for mid, met in new_mets.items():
        if mid not in existing:
            model.add_metabolites([met])
            print(f"  Added metabolite: {mid}")
        else:
            print(f"  Metabolite already present: {mid}")
    
    # ── transporter: hdec_e → hdec_c ─────────────────────────────────────────
    transport = Reaction("HDECtex")
    transport.name        = "Hexadecanoic acid transporter (extracellular to cytosol)"
    transport.lower_bound = 0
    transport.upper_bound = 1000
    transport.add_metabolites({
        model.metabolites.get_by_id("hdec_e"): -1.0,
        model.metabolites.get_by_id("hdec_c"):  1.0,
    })
    model.add_reactions([transport])
    
    # ── exchange reaction for hdec_e ──────────────────────────────────────────
    ex_hdec = Reaction("EX_hdec_e")
    ex_hdec.name        = "Hexadecanoic acid exchange"
    ex_hdec.lower_bound = -1000
    ex_hdec.upper_bound =  1000
    ex_hdec.add_metabolites({model.metabolites.get_by_id("hdec_e"): -1.0})
    model.add_reactions([ex_hdec])
    
    # ── palmitate CoA ligase: atp + coa + h + hdec → amp + ppi + pmtcoa ──────
    palmcl = Reaction("PALMCL")
    palmcl.name        = "Palmitate CoA ligase"
    palmcl.lower_bound = 0
    palmcl.upper_bound = 1000
    palmcl.add_metabolites({
        model.metabolites.get_by_id("atp_c")    : -1.0,
        model.metabolites.get_by_id("coa_c")    : -1.0,
        model.metabolites.get_by_id("h_c")      : -1.0,
        model.metabolites.get_by_id("hdec_c")   : -1.0,
        model.metabolites.get_by_id("amp_c")    :  1.0,
        model.metabolites.get_by_id("ppi_c")    :  1.0,
        model.metabolites.get_by_id("pmtcoa_c") :  1.0,
    })
    model.add_reactions([palmcl])
    
    # ── remove peroxisome-only reactions ─────────────────────────────────────
    to_remove = ["PGCD", "PSERT", "CS", "DPCOAK", "SACCD2"]
    model.remove_reactions([model.reactions.get_by_id(r) for r in to_remove])
    print(f"  Removed peroxisome reactions: {to_remove}")
    
def modify_original_model_to_updated_no_peroxisomes(model):
    """All structural changes specific to the updated model."""
    
    # ── check which metabolites already exist ─────────────────────────────────
    existing = [m.id for m in model.metabolites]
    new_mets = {
        "hdec_e"  : Metabolite("hdec_e",   name="Hexadecanoic acid (extracellular)", compartment="e"),
        "hdec_c"  : Metabolite("hdec_c",   name="Hexadecanoic acid",                compartment="c"),
        "pmtcoa_c": Metabolite("pmtcoa_c", name="Palmitoyl-CoA",                    compartment="c"),
        "ppi_c"   : Metabolite("ppi_c",    name="Diphosphate",                      compartment="c"),
        "amp_c"   : Metabolite("amp_c",    name="AMP",                              compartment="c"),
    }
    for mid, met in new_mets.items():
        if mid not in existing:
            model.add_metabolites([met])
            print(f"  Added metabolite: {mid}")
        else:
            print(f"  Metabolite already present: {mid}")
    
    # ── transporter: hdec_e → hdec_c ─────────────────────────────────────────
    transport = Reaction("HDECtex")
    transport.name        = "Hexadecanoic acid transporter (extracellular to cytosol)"
    transport.lower_bound = 0
    transport.upper_bound = 1000
    transport.add_metabolites({
        model.metabolites.get_by_id("hdec_e"): -1.0,
        model.metabolites.get_by_id("hdec_c"):  1.0,
    })
    model.add_reactions([transport])
    
    # ── exchange reaction for hdec_e ──────────────────────────────────────────
    ex_hdec = Reaction("EX_hdec_e")
    ex_hdec.name        = "Hexadecanoic acid exchange"
    ex_hdec.lower_bound = -1000
    ex_hdec.upper_bound =  1000
    ex_hdec.add_metabolites({model.metabolites.get_by_id("hdec_e"): -1.0})
    model.add_reactions([ex_hdec])
    
    # ── palmitate CoA ligase: atp + coa + h + hdec → amp + ppi + pmtcoa ──────
    palmcl = Reaction("PALMCL")
    palmcl.name        = "Palmitate CoA ligase"
    palmcl.lower_bound = 0
    palmcl.upper_bound = 1000
    palmcl.add_metabolites({
        model.metabolites.get_by_id("atp_c")    : -1.0,
        model.metabolites.get_by_id("coa_c")    : -1.0,
        model.metabolites.get_by_id("h_c")      : -1.0,
        model.metabolites.get_by_id("hdec_c")   : -1.0,
        model.metabolites.get_by_id("amp_c")    :  1.0,
        model.metabolites.get_by_id("ppi_c")    :  1.0,
        model.metabolites.get_by_id("pmtcoa_c") :  1.0,
    })
    model.add_reactions([palmcl])
    
    # ── beta oxidation ──────
    FAOX = Reaction("FAOX")
    FAOX.name        = "Beta oxidation"
    FAOX.lower_bound = 0
    FAOX.upper_bound = 1000
    FAOX.add_metabolites({
        model.metabolites.get_by_id("coa_c")       : -7.0,
        model.metabolites.get_by_id("fad_c")       : -7.0,
        model.metabolites.get_by_id("h2o_c")       : -7.0,
        model.metabolites.get_by_id("nad_c")       : -7.0,
        model.metabolites.get_by_id("pmtcoa_c")    : -1.0,
        model.metabolites.get_by_id("accoa_c")     :  8.0,
        model.metabolites.get_by_id("fadh2_c")     :  7.0,
        model.metabolites.get_by_id("nadh_c")      :  7.0,
        model.metabolites.get_by_id("h_c")         :  7.0,
    })
    model.add_reactions([FAOX])
    
    # ── OSM1 fad dependent ──────
    OSM = Reaction("OSM")
    OSM.name        = "Beta oxidation"
    OSM.lower_bound = 0
    OSM.upper_bound = 1000
    OSM.add_metabolites({
        model.metabolites.get_by_id("fum_c")         : -1.0,
        model.metabolites.get_by_id("fadh2_c")       : -1.0,
        model.metabolites.get_by_id("succ_c")        :  1.0,
        model.metabolites.get_by_id("fad_c")         :  1.0,
    })
    model.add_reactions([OSM])
    
def activate_selected_pathways(model, pathways):
    """Turn on selected pathways."""
    for rxn_id in pathways:
        rxn             = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = 0
        rxn.upper_bound = 1000
    print("  Activated selected pathways")
    
def deactivate_selected_pathways(model, pathways):
    """Turn off selected pathways."""
    for rxn_id in pathways:
        rxn             = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    print("  Deactivated selected pathways")
    
############################
## function to check FAOX ##
############################
def run_reduce_stepswise(model_key, rxn_id, runSample = False, update_original = False, update = True, constrains = True, steps = [0.005, 0.004, 0.003, 0.002, 0.001,  0]):
    """Runs my updated model and slowly turns of the faox and returs fluxes"""
    #load the model
    print(f"  RUNNING: {model_key} model")
    model = cobra.io.read_sbml_model(MODEL_PATHS[model_key])
    model.solver = "gurobi"
    
    #prepare the model
    if update:
        print("Updating the model")
        modify_updated_model(model)
        activate_selected_pathways(model, ["HYDh2", "HYDh3"]) # activation hydrogenases
    
    if update_original:
        print("Updating the original model to mimic the updated without peroxisomes")
        modify_original_model_to_updated_no_peroxisomes(model)
        activate_selected_pathways(model, ["HYDh2", "HYDh3"]) # activation hydrogenases
    
    if constrains:
        apply_exchange_constraints(model, EXCHANGE_CONSTRAINTS) #applay constrains
    
    #define steps of reduction
    print(f"Reducing the flux of {rxn_id} and runing the ")
    
    for step in steps:
        #make the tag
        step_ch = str(step)
        tag_nocon = f"{model_key}_{rxn_id}_{step_ch}"
        
        #reduce the flux
        rxn             = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = 0
        rxn.upper_bound = step
        
        #run the model
        run_fba(model, tag_nocon)
        
        if runSample:           
            run_sampling(model, tag_nocon)
    
####################
## functions that go through the constrains and see which combination of constrains is fiseable 
###################
def check_constraint_combinations(model, faox_flux=0.001,
                                   constraints=EXCHANGE_CONSTRAINTS,
                                   flux_reactions=["PFOh", "PFLh", "LDH_D", "HYDhfe",
                                   "HYDh3", "HYDh2", "SUCOAACTh", "ALCD2y", "ALCD2x",
                                   "ME1", "ME2", "MDHh", "AMTh", "MAKht", "AGht",
                                   "FADxt", "COAxt", "MOAAxt", "AKGxt", "ASPxt",
                                   "FAOXx", "ASTx", "MDHx", "FUMh", "PYRht",
                                   "EX_h2_e", "EX_ac_e", "EX_for_e", "EX_h_e", "EX_co2_e",
                                   "EX_etoh_e", "EX_lac__D_e", "Biomass", "EX_succ_e",
                                   "ALDD2x", "FRDx",
                                   ],
                                   verbose=True):
    rxn_ids = list(constraints.keys())
    modes   = ["on", "off"]
    results = []
    
    combos = list(product(modes, repeat=len(rxn_ids)))
    print(f"Testing {len(combos)} combinations...")
    
    for combo in combos:
        config = dict(zip(rxn_ids, combo))
        applied_bounds = {}
        
        with model:
            faox = model.reactions.get_by_id("FAOXx")
            faox.lower_bound = faox_flux
            faox.upper_bound = faox_flux
            
            for rxn_id, mode in config.items():
                rxn = model.reactions.get_by_id(rxn_id)
                if mode == "on":
                    lb, ub = constraints[rxn_id]
                    rxn.lower_bound = lb
                    rxn.upper_bound = ub
                    applied_bounds[rxn_id] = (lb, ub)
                else:  # "off"
                    rxn.lower_bound = -1000
                    rxn.upper_bound = 1000
                    applied_bounds[rxn_id] = (-1000, 1000)
                
            sol = model.optimize()
            status = sol.status
        
        row = {
            "config" : config,
            "bounds" : applied_bounds,
            "status" : status,
        }
        
        # Add fluxes only when feasible
        if status == "optimal":
            for rxn_id in flux_reactions:
                try:
                    row[rxn_id] = sol.fluxes[rxn_id]
                except KeyError:
                    row[rxn_id] = None   # reaction not in model
        else:
            for rxn_id in flux_reactions:
                row[rxn_id] = float("nan")
        
        results.append(row)
        
        if verbose:
            feasible = "✓" if status == "optimal" else "✗"
            cfg_str  = " | ".join(
                f"{k.replace('EX_','').replace('_e','')}:{v}"
                for k, v in config.items()
            )
            print(f"  {feasible}  {cfg_str}  →  {status}")
    
    df = pd.DataFrame(results)
    n_feasible   = (df["status"] == "optimal").sum()
    n_infeasible = (df["status"] == "infeasible").sum()
    print(f"\nSummary: {n_feasible} feasible / {n_infeasible} infeasible "
          f"out of {len(combos)} combinations")
    return df
    
##############################
## main pipeline            ##
##############################
def run_pipeline(model_key):
    print(f"\n{'='*60}")
    print(f"  RUNNING: {model_key} model")
    print(f"{'='*60}")
    
    model = cobra.io.read_sbml_model(MODEL_PATHS[model_key])
    model.solver = "gurobi"
    
    #use cellebiose instead of glucose
    glc_ex = model.reactions.get_by_id("EX_glc__D_e")
    glc_ex.upper_bound = 0
    glc_ex.lower_bound = 0
    
    #pripadne zmencit intake na pul since pu hydrolasa vyrobi 2 glucosy z 1 cellebiosy
    cel_ex = model.reactions.get_by_id("EX_cellb_e")
    cel_ex.lower_bound = -0.75
    cel_ex.upper_bound = -0.75
    
    # model-specific modifications
    if model_key == "updated":
        modify_updated_model(model)
        activate_selected_pathways(model, ["HYDh2", "HYDh3"]) # activation hydrogenases
    elif model_key == "updated_noFAOX":
        modify_updated_model(model)
        activate_selected_pathways(model, ["HYDh2", "HYDh3"]) # activation hydrogenases
        deactivate_selected_pathways(model, ["FAOXx"]) # deactivate FAOX
    
    # ── no-constraints run ────────────────────────────────────────────────────
    print("\n-- No constraints --")
    tag_nocon = f"{model_key}_NoCon"
    run_fba(model, tag_nocon)
    run_fva(model, tag_nocon)
    run_sampling(model, tag_nocon)
    save_model_tables(model, model_key)
    
    # ── with constraints run ──────────────────────────────────────────────────
    print("\n-- With constraints --")
    apply_exchange_constraints(model, EXCHANGE_CONSTRAINTS)
    tag_withcon = f"{model_key}_WithCon"
    fba_sol = run_fba(model, tag_withcon)
    run_fva(model, tag_withcon)
    apply_biomass_constrain(model, fba_sol)
    run_sampling(model, tag_withcon)

####################
## run the models ##
####################
if __name__ == "__main__":
    run_pipeline("original")
    run_pipeline("updated")
#    run_pipeline("updated_noFAOX") #cannot pass, it fails
    
    #####################
    ## limit FAOX flux ##
    #####################
    #shows the decrease in biomass
    run_reduce_stepswise("updated", rxn_id = "FAOXx", update_original = False, update = True, constrains = True, steps = [0.005, 0.004, 0.003, 0.002, 0.001,  0]) 
