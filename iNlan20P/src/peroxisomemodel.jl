#################################################################
# Model of the peroxisome
#################################################################

##########################################
# Metabolites
##########################################
# must be in BIGG-style standardized metabolite nomenclature.
# the _x means peroxisomal

#voda
met = GSM.createmetabolite("h2o_x", unimets)
GSM.addmet!(model, met)

#dephospho coa
met = GSM.createmetabolite("dpcoa_x", unimets)
GSM.addmet!(model, met)

#coa
met = GSM.createmetabolite("coa_x", unimets)
GSM.addmet!(model, met)

#acetyl-coa
met = GSM.createmetabolite("accoa_x", unimets)
GSM.addmet!(model, met)

#atp
met = GSM.createmetabolite("atp_x", unimets)
GSM.addmet!(model, met)

#adp
met = GSM.createmetabolite("adp_x", unimets)
GSM.addmet!(model, met)

#gtp
met = GSM.createmetabolite("gtp_x", unimets)
GSM.addmet!(model, met)

#gdp
met = GSM.createmetabolite("gdp_x", unimets)
GSM.addmet!(model, met)

#co2
met = GSM.createmetabolite("co2_x", unimets)
GSM.addmet!(model, met)

#nad
met = GSM.createmetabolite("nad_x", unimets)
GSM.addmet!(model, met)

#nadh
met = GSM.createmetabolite("nadh_x", unimets)
GSM.addmet!(model, met)

#fad
met = GSM.createmetabolite("fad_x", unimets)
GSM.addmet!(model, met)

#fadh2
met = GSM.createmetabolite("fadh2_x", unimets)
GSM.addmet!(model, met)

#fummarate
met = GSM.createmetabolite("fum_x", unimets)
GSM.addmet!(model, met)

#succinate
met = GSM.createmetabolite("succ_x", unimets)
GSM.addmet!(model, met)

#l-lysine
met = GSM.createmetabolite("lys__L_x", unimets)
GSM.addmet!(model, met)

#oaa
met = GSM.createmetabolite("oaa_x", unimets)
GSM.addmet!(model, met)

#pep
met = GSM.createmetabolite("pep_x", unimets)
GSM.addmet!(model, met)

#citrate
met = GSM.createmetabolite("cit_x", unimets)
GSM.addmet!(model, met)

#l-glutamate
met = GSM.createmetabolite("glu__L_x", unimets)
GSM.addmet!(model, met)

#l-asparate
met = GSM.createmetabolite("asp__L_x", unimets)
GSM.addmet!(model, met)

#alpha-ketoglutarate
met = GSM.createmetabolite("akg_x", unimets)
GSM.addmet!(model, met)

#l-malate
met = GSM.createmetabolite("mal__L_x", unimets)
GSM.addmet!(model, met)

#3-Phosphohydroxypyruvate
met = GSM.createmetabolite("3php_x", unimets)
GSM.addmet!(model, met)

#3-Phospho-D-glycerate
met = GSM.createmetabolite("3pg_x", unimets)
GSM.addmet!(model, met)

#3P-L-serine
met = GSM.createmetabolite("pser__L_x", unimets)
GSM.addmet!(model, met)

#L Saccharopine 
met = GSM.createmetabolite("saccrp__L_x", unimets)
GSM.addmet!(model, met)

##D-Fructose 6-phosphate
#met = GSM.createmetabolite("f6p_x", unimets)
#GSM.addmet!(model, met)

##D-Fructose 1,6-bisphosphate
#met = GSM.createmetabolite("fdp_x", unimets)
#GSM.addmet!(model, met)

##glyceraldehyd-3p
#met = GSM.createmetabolite("g3p_x", unimets)
#GSM.addmet!(model, met)

##Alpha-D-Ribose 5-phosphate
#met = GSM.createmetabolite("r5p_x", unimets)
#GSM.addmet!(model, met)

#hydrogen proton
met = GSM.createmetabolite("h_x", unimets)
GSM.addmet!(model, met)

##D-Xylulose 5-phosphate
#met = GSM.createmetabolite("xu5p__D_x", unimets)
#GSM.addmet!(model, met)

##Sedoheptulose 7-phosphate
#met = GSM.createmetabolite("s7p_x", unimets)
#GSM.addmet!(model, met)

##D-Erythrose 4-phosphate
#met = GSM.createmetabolite("e4p_x", unimets)
#GSM.addmet!(model, met)

#######################
##betaoxidation je brana jako 1 reakce ti padem nepotrebuju meziprodukty a davam pouze zacatecni fatty acid a pak tu o 2 C kratsi
#######################
#mozna musim dodat ale vsechny produkty az po posledni C4 acyl!!!!!!!!

#Hexadecenoyl-CoA (n-C16:1CoA)
met = GSM.createmetabolite("pmtcoa_x", unimets)
GSM.addmet!(model, met)

#################################################################
# Transporters
#################################################################
####################
## unidirectional ##
####################
# dephospho coa transporter
md = Dict("dpcoa_c" => -1, "dpcoa_x" => 1)
rxn = GSM.createreaction("DPCOAxt", "3'-dephospho-CoA peroxisome transporter", md)
GSM.addrxn!(model, rxn)

# Fatty acyl-CoA import (ABCD transporter) (Palmitoyl-CoA)
md = Dict("pmtcoa_c" => -1, "pmtcoa_x" => 1)
rxn = GSM.createreaction("ABCDxt", "Hexadecenoyl-CoA peroxisome ABCD transporter", md)
GSM.addrxn!(model, rxn)

# 3P-L-serine transport
md = Dict("pser__L_x" => -1, "pser__L_c" => 1)
rxn = GSM.createreaction("PSERLxt", "3P-L-serine peroxisome transporter", md)
GSM.addrxn!(model, rxn)

###################
## bidirectional ##
###################
# diffusion h2o
md = Dict("h2o_c" => -1, "h2o_x" => 1)
rxn = GSM.createreaction("H2Oxt", "H2O peroxisome transporter (diffusion)", md)
GSM.addrxn!(model, rxn)

# diffusion co2
md = Dict("co2_c" => -1, "co2_x" => 1)
rxn = GSM.createreaction("CO2xt", "CO2 peroxisome transporter (diffusion)", md)
GSM.addrxn!(model, rxn)

# diffusion h+
md = Dict("h_c" => -1, "h_x" => 1)
rxn = GSM.createreaction("Hxt", "Proton peroxisome transporter (diffusion)", md)
GSM.addrxn!(model, rxn)

# ATP/ADP antiport
md = Dict("atp_c" => -1, "atp_x" => 1, "adp_c" => 1, "adp_x" => -1)
rxn = GSM.createreaction("ATPADPxt", "ATP/ADP peroxisome antiporter", md)
GSM.addrxn!(model, rxn)

# GTP/GDP antiport
md = Dict("gtp_c" => -1, "gtp_x" => 1, "gdp_c" => 1, "gdp_x" => -1)
rxn = GSM.createreaction("GTPGDPxt", "GTP/GDP peroxisome antiporter", md)
GSM.addrxn!(model, rxn)

#not needed that is done by maltase shuttle
## NAD+/NADH antiport = netreba mel by to zvladat OAA malate shuttle
#md = Dict("nad_c" => -1, "nad_x" => 1, "nadh_c" => 1, "nadh_x" => -1)
#rxn = GSM.createreaction("NADxt", "NAD+/NADH peroxisome antiporter", md)
#GSM.addrxn!(model, rxn)

#not true
## FAD/FADH2 antiport
#md = Dict("fad_c" => -1, "fad_x" => 1, "fadh2_c" => 1, "fadh2_x" => -1)
#rxn = GSM.createreaction("FADxt", "FAD/FADH2 peroxisome antiporter", md)
#GSM.addrxn!(model, rxn)

#possibly by pmp34
# FAD transport namisto toho FAD/FADH2 antiport
md = Dict("fad_c" => -1, "fad_x" => 1)
rxn = GSM.createreaction("FADxt", "FAD peroxisome transport", md)
GSM.addrxn!(model, rxn)

#possibly by pmp34
# CoA transporter
md = Dict("coa_c" => -1, "coa_x" => 1)
rxn = GSM.createreaction("COAxt", "CoA peroxisome transporter", md)
GSM.addrxn!(model, rxn)

## not possible by pmp34?
#md = Dict("accoa_c" => -1, "accoa_x" => 1)
#rxn = GSM.createreaction("ACCOAxt", "Acetyl-CoA peroxisome transporter", md)
#GSM.addrxn!(model, rxn)

# Malate/OAA antitransport
md = Dict("mal__L_c" => 1, "mal__L_x" => -1, "oaa_c" => -1, "oaa_x" => 1)
rxn = GSM.createreaction("MOAAxt", "Malate/OAA peroxisome antitransporter", md)
GSM.addrxn!(model, rxn)

# Alpha-ketoglutarate transport
md = Dict("akg_c" => -1, "akg_x" => 1)
rxn = GSM.createreaction("AKGxt", "Alpha-ketoglutarate peroxisome transporter", md)
GSM.addrxn!(model, rxn)

# Citrate transport
md = Dict("cit_c" => 1, "cit_x" => -1)
rxn = GSM.createreaction("CITxt", "Citrate peroxisome transporter", md)
GSM.addrxn!(model, rxn)

# 3-Phospho-D-glycerate transport
md = Dict("3pg_c" => -1, "3pg_x" => 1)
rxn = GSM.createreaction("3PGxt", "3-Phospho-D-glycerate peroxisome transporter", md)
GSM.addrxn!(model, rxn)

# PEP
md = Dict("pep_c" => -1, "pep_x" => 1)
rxn = GSM.createreaction("PEPxt", "Phosphoenolpyruvate peroxisome transporter", md)
GSM.addrxn!(model, rxn)

# Succinate transport
md = Dict("succ_c" => -1, "succ_x" => 1)
rxn = GSM.createreaction("SUCCxt", "Succinate peroxisome transporter", md)
GSM.addrxn!(model, rxn)

# Fumarate transport
md = Dict("fum_c" => -1, "fum_x" => 1)
rxn = GSM.createreaction("FUMxt", "Fumarate peroxisome transporter", md)
GSM.addrxn!(model, rxn)

# Saccharopine transport
md = Dict("saccrp__L_c" => -1, "saccrp__L_x" => 1)
rxn = GSM.createreaction("SACCRPxt", "Saccharopine peroxisome transporter", md)
GSM.addrxn!(model, rxn)

# Lysine transport
md = Dict("lys__L_c" => -1, "lys__L_x" => 1)
rxn = GSM.createreaction("LYSLxt", "Lysine peroxisome transporter", md)
GSM.addrxn!(model, rxn)

## ATP/ADP symport
#md = Dict("atp_c" => -1, "atp_x" => 1, "adp_c" => -1, "adp_x" => 1)
#rxn = GSM.createreaction("ATPxt", "ATP/ADP peroxisome symporter", md)
#GSM.addrxn!(model, rxn)

## GTP/GDP symport???????
#md = Dict("gtp_c" => -1, "gtp_x" => 1, "gdp_c" => -1, "gdp_x" => 1)
#rxn = GSM.createreaction("GTPxt", "GTP/GDP peroxisome symporter", md)
#GSM.addrxn!(model, rxn)

# Aspartate glutamate antitransport
md = Dict("asp__L_c" => -1, "asp__L_x" => 1, "glu__L_c" => 1, "glu__L_x" => -1)
rxn = GSM.createreaction("ASPxt", "Aspartate/glutamate peroxisome antitransporter", md)
GSM.addrxn!(model, rxn)


######################################################
## reactions 
######################################################
###################################
## unidirectional 
###################################

###################
# Beta oxidation ##
###################
md = Dict("coa_x" => -7, 
		  "fad_x" => -7, 
		  "h2o_x" => -7,
		  "nad_x" => -7,
		  "pmtcoa_x" => -1,
		  "accoa_x" => 8,
		  "fadh2_x" => 7,
		  "nadh_x" => 7,
		  "h_x" => 7
		  )
rxn = GSM.createreaction("FAOXx", "Beta oxidation peroxisome", md)
GSM.addrxn!(model, rxn)

##########
## DPCK ##
##########
md = Dict("atp_x" => -1,
		  "dpcoa_x" => -1,
		  "adp_x" => 1,
		  "coa_x" => 1,
		  "h_x" => 1
		  )
rxn = GSM.createreaction("DPCKx", "Generation of CoA for beta-oxidation peroxisome", md)
GSM.addrxn!(model, rxn)

######################
## citrate synthase ##
######################
md = Dict("oaa_x" => -1,
		  "accoa_x" => -1,
		  "h2o_x" => -1,
		  "cit_x" => 1,
		  "coa_x" => 1,
		  "h_x" => 1
		  )
rxn = GSM.createreaction("CSx", "Citrate synthase peroxisome", md)
GSM.addrxn!(model, rxn)

######################################
## succinate to fumarate conversion ##
######################################
md = Dict("fum_x" => -1,
		  "fadh2_x" => -1,
		  "succ_x" => 1,
		  "fad_x" => 1
		  )
rxn = GSM.createreaction("OSMx", "Succinate to Fumarate conversion peroxisome", md)
GSM.addrxn!(model, rxn)

##########
## PGDH ##
##########
md = Dict("3pg_x" => -1,
		  "nad_x" => -1,
		  "3php_x" => 1,
		  "nadh_x" => 1,
		  "h_x" => -1
		  )
rxn = GSM.createreaction("PGDHx", "Phosphoglycerate dehydrogenase peroxisome", md)
GSM.addrxn!(model, rxn)

###############################
## bidirectional 
###############################
##########################
## oaa synthesis of aKG ##
##########################
md = Dict("oaa_x" => 1,
		  "glu__L_x" => 1,
		  "akg_x" => -1,
		  "asp__L_x" => -1
		  )
rxn = GSM.createreaction("ASTx", "Aspartate aminotransferase peroxisome", md)
GSM.addrxn!(model, rxn)

#####################################
## MDH coupled with beta oxidation ##
#####################################
md = Dict("oaa_x" => -1,
		  "nadh_x" => -1,
		  "h_x" => -1,
		  "nad_x" => 1,
		  "mal__L_x" => 1
		  )
rxn = GSM.createreaction("MDHx", "Malate Dehydrogenase peroxisome", md)
GSM.addrxn!(model, rxn)

#####################
## saccharopine DH ##
#####################
md = Dict("saccrp__L_x" => -1,
		  "nad_x" => -1,
		  "h2o_x" => -1,
		  "akg_x" => 1,
		  "nadh_x" => 1,
		  "h_x" => 1,
		  "lys__L_x" => 1
		  )
rxn = GSM.createreaction("SADHx", "Succinate to Fumarate conversion peroxisome", md)
GSM.addrxn!(model, rxn)

###########
## PEPCK ##
###########
md = Dict("oaa_x" => -1,
		  "gdp_x" => 1,
		  "co2_x" => 1,
		  "gtp_x" => -1,
		  "pep_x" => 1
		  )
rxn = GSM.createreaction("PEPCKx", "Phosphoenolpyruvate carboxykinase peroxisome", md)
GSM.addrxn!(model, rxn)

##########
## PSAT ##
##########
md = Dict("3php_x" => -1,
		  "glu__L_x" => -1,
		  "akg_x" => 1,
		  "pser__L_x" => 1
		  )
rxn = GSM.createreaction("PSATx", "Phosphoserine aminotransferase peroxisome", md)
GSM.addrxn!(model, rxn)

#########
### TK ##
#########
#md = Dict("r5p_x" => -1,
#		  "xu5p__D_x" => -1,
#		  "g3p_x" => 1,
#		  "s7p_x" => 1
#		  )
#rxn = GSM.createreaction("TKx", "Transketolase peroxisome", md)
#GSM.addrxn!(model, rxn)
#
#########
### TA ##
#########
#md = Dict("g3p_x" => -1,
#		  "s7p_x" => -1,
#		  "f6p_x" => 1,
#		  "e4p_x" => 1
#		  )
#rxn = GSM.createreaction("TAx", "Transaldolase peroxisome", md)
#GSM.addrxn!(model, rxn)
#
##########
### PFK ##
##########
#md = Dict("f6p_x" => -1,
#		  "atp_x" => -1,
#		  "adp_x" => 1,
#		  "fdp_x" => 1
#		  )
#rxn = GSM.createreaction("PFKx", "Phosphofructokinase  peroxisome", md)
#GSM.addrxn!(model, rxn)

