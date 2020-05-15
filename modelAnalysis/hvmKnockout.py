import  os
from    os.path     import join
import  cobra
from    cobra       import Model, Reaction, Metabolite
from cobra.flux_analysis.loopless import loopless_solution
from cobra.util.solver import linear_reaction_coefficients
import  numpy       as np
import  scipy       as sp
import  csv         as cp
import pandas as pd
import  sys
import  re
#######################################################################################
# Function Definition
# hvmKnockout Evaluates the effect of enforcing a knockout state on virus
# production (model generated by genHVM.py)

# Inputs:
# HVM               Integrated host-virus model
# HostRxn           Name of the host objective reaction

# Outputs:
# koVirus          Vector of virus optima values with additional host-constraint
def hvmKnockout(HVM,HostRxn,solver):
	"""Reaction Knockouts
	new version by Hadrien"""

	# [1] Optimise the HVM for host and virus objectives
	# Initial Setup
	HVM.solver = solver										#set solver for all optimisations
	hostIdx = next(index for index, reaction in enumerate(HVM.reactions) if reaction.id == HostRxn)
	virusIdx    = len(HVM.reactions) - 1
	objIdx      = [hostIdx, virusIdx]
	virusRxn = HVM.reactions[-1].id

	# Host Optimisation
	hostObj     = HVM.reactions[hostIdx]
	HVM.objective   = hostObj.id
	hostSol     = HVM.optimize()
	print("Host optimization objective: {}".format(hostSol.objective_value))

	# Virus Optimisation
	virusObj    = HVM.reactions[-1]
	HVM.objective   = virusObj.id
	virusSol    = HVM.optimize()
	print("Virus  optimization objective: {}".format(virusSol.objective_value))

	# [2] Knockout Analysis
	koVirus             = np.zeros((len(HVM.reactions),5))

	# Initiate loop
	nbofReactions = len(HVM.reactions)
	reactionNames = [None] * nbofReactions
	for ii in range(nbofReactions):
		# Conditional to exclude objective reactions and those that carry zero flux for virus optima
		if ((ii!=hostIdx) and (ii!=virusIdx) and (virusSol.fluxes[ii]!=0.0) and (not HVM.reactions[ii].id.startswith("EX_"))):
			print("{}/{}\tKnockout: {}\t".format(ii, nbofReactions, HVM.reactions[ii].id))
			# Store reaction info
			reactionNames[ii] = HVM.reactions[ii].id
			koVirus[ii,0] = ii
			# Store original fluxes under host and virus optimisation
			koVirus[ii,1] = hostSol.fluxes[ii]
			koVirus[ii,2] = virusSol.fluxes[ii]
			# Store the bounds
			tempLB = HVM.reactions[ii].lower_bound
			tempUB = HVM.reactions[ii].upper_bound
			# Alter bounds to zero
			HVM.reactions[ii].lower_bound = 0
			HVM.reactions[ii].upper_bound = 0
			# Optimise the model for virus or host production
			HVM.objective   = virusObj.id  
			koVirus[ii,3] = HVM.slim_optimize()
			HVM.objective   = hostObj.id    
			koVirus[ii,4] = HVM.slim_optimize()
			# Return to original bounds
			HVM.reactions[ii].lower_bound = tempLB
			HVM.reactions[ii].upper_bound = tempUB
		else:
			koVirus[ii,3] = np.nan
			koVirus[ii,4] = np.nan
	print("done")

	# [3] Post Analysis
	# Remove nan values
	# Convert to % original optima
	for ii in range(len(koVirus)):
		koVirus[ii,3] = (koVirus[ii,3] / virusSol.objective_value) * 100
		koVirus[ii,4] = (koVirus[ii,4] / hostSol.objective_value) * 100

	# [4] Output
	outputDf = pd.DataFrame(data=koVirus, columns=["index", "host_flux", "virus_flux", "virus_optima_KO","host_optima_KO"])
	outputDf = outputDf.assign(name=reactionNames)
	return outputDf
