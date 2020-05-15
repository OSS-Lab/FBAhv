########################### Identification of a Minimal Medium  ###########################
###########################################################################################
# Written by: Kalesh Sasidharan
# Email: k.sasidharan@warwick.ac.uk
# :: Released on: 2020 April 08
###########################################################################################

# Setting up the CobraPy
from __future__ import print_function
import cobra
from cobra import Reaction, Metabolite
import numpy as np
import scipy

###########################################################################################
# A function for finding a Minimal Medium
# Required functions: showExcRxns, restrictBounds, addC, findHighYieldImport, findCritical, uptakeEliminator, findActiveUptake
# Algorithm:
# 1. Find high-yield uptake reactions using an algorithm similar to the Membrane Economics algorithm. See "findHighYieldImport" for details.
# 2. Find the critical uptake reactions among the high-yield uptake reactions.
#    If the critical uptake reactions can produce a desirable growth rate (provided via the "desirableGrowth" parameter), then that is the best minimal medium.
# 3. Eliminate the non-critical high-yield uptake reactions that have minimal effect on the growth rate.
# 4. Check the Minimal Medium uptake reactions to make sure that they all carry flux. If not, select the ones carrying flux.
# 5. If possible, it returns a Minimal Medium that meets the minimum growth rate and "omit" requirements.
# See "findHighYieldImport" for details about "model", "omit" and "growthCutOff" parameters.
# The "desirableGrowth" parameter determines the minimum growth rate that should be produced by a Minimal Medium.
# Return values: indexes of the Minimal Medium exchange reactions in the model, indexes of critical Minimal Medium reactions, flux value of Minimal Medium reactions,
# and optimisation result using Minimal Medium.
def findMinimalMedium(model, omit, growthCutOff, desirableGrowth):
  rxn = []
  rxnInd =[]
  omitInd = []
  hyInd = []
  activeUptake = []
  tmpModel = model.copy()
  # Finding the uptake reactions
  modelExcUps = showExcRxns(tmpModel)
  # Finding the high-yield uptake reactions
  print("\nFinding the high-yield uptake reactions.\n")
  highYieldImport = findHighYieldImport(tmpModel, omit, growthCutOff)
  # Resetting the model
  tmpModel = model.copy()
  # Finding the critical uptake reactions
  print("\nFinding the critical uptake reactions.\n")
  criticalUptake = findCritical(tmpModel, highYieldImport[0], growthCutOff)
  # Restricting all the uptake reactions
  restrictBounds(tmpModel, modelExcUps[3], 0, 1000)
  # Applying only the critical uptake reactions
  restrictBounds(tmpModel, criticalUptake[0], -1000, 1000)
  # If the critical components can produce the desirable growth then that is the best minimal medium
  if (tmpModel.slim_optimize() >= desirableGrowth):
    print ("\nAs the model can be optimised using only the critical media components, returning those exchange reactions. Note that this list may not include the 'omit' reactions, if you passed any.\n")
    return criticalUptake[0]
  # Resetting the model
  tmpModel = model.copy()
  # Restricting all the uptake reactions
  restrictBounds(tmpModel, modelExcUps[3], 0, 1000)
  # Applying the high-yield medium
  for rxnId in highYieldImport[0]:
    rxn = tmpModel.reactions.get_by_id(rxnId)
    rxnInd.append(tmpModel.reactions.index(rxn))
  restrictBounds(tmpModel, rxnInd, -1000, 1000)
  # Finding the high-yield medium components that are not critical
  if (omit != ""):
    for rxnId in omit:
      rxn = tmpModel.reactions.get_by_id(rxnId)
      omitInd.append(tmpModel.reactions.index(rxn))
    omitInd = list(set(rxnInd).intersection(set(omitInd)))
  hyInd = list(set(rxnInd) - set(criticalUptake[0]) - set(omitInd))
  # Eliminating the unnecessary uptake reactions
  print("\nEliminating the unnecessary uptake reactions.")
  minimalMedia = uptakeEliminator(tmpModel, hyInd, desirableGrowth)
  minimalMedia = minimalMedia + criticalUptake[0] + omitInd
  # Checking the Minimal Medium uptake reactions to make sure that they all are still carrying flux. If not, select the ones carrying flux.
  print("\nTesting and optimising the Minimal Medium.\n")
  # Resetting the model
  tmpModel = model.copy()
  activeUptake = findActiveUptake(tmpModel, minimalMedia, growthCutOff, 1)
  return activeUptake[0], criticalUptake[0], activeUptake[1], activeUptake[2]
###########################################################################################

###########################################################################################
# A function for identifying the exchange and uptake reactions, their IDs and indices
def showExcRxns(model):
  exIndex = []
  exID = []
  exNames = []
  upIndex = []
  upID = []
  upNames = []
  upMetID = []
  upMetNames = []
  for rxn in model.reactions:
    # Exchange reactions
    if rxn.products == [] or rxn.reactants == []:
      exIndex.append(model.reactions.index(rxn))
      exID.append(rxn.id)
      exNames.append(rxn.name)
      # Uptake reactions
      if rxn.lower_bound < 0:
        upIndex.append(model.reactions.index(rxn))
        upID.append(rxn.id)
        upNames.append(rxn.name)
        upMetID.append(list(rxn.metabolites)[0].id)
        upMetNames.append(list(rxn.metabolites)[0].name)
  return exIndex, exID, exNames, upIndex, upID, upNames, upMetID, upMetNames
###########################################################################################

###########################################################################################
# Restrict all the given reactions
def restrictBounds(model, indices, lb, ub):
  for ind in indices:
    rxn = model.reactions[ind]
    rxn.upper_bound = ub
    rxn.lower_bound = lb
###########################################################################################

###########################################################################################
# A function for adding an imaginary metabolite "c" to all the uptake reactions (e.g., x[e] <=>  to x[e] + c <=>).
# Export reactions are excluded. Therefore, bidirectional exchange reactions are split into uptake and export reactions.
# One can omit reaction(s) from adding 'c' by passing the reaction IDs as an array via the "omit" parameter.
def addC(model, omit):
  fluxLimit = 0
  rxnTMP = []
  # Finding the exchange/uptake reactions
  modelExcUps = showExcRxns(model)
  # Finding the total flux carried by all uptake reactions for setting initial upper bound for the "c" export reaction
  modelFlux = model.optimize()
  for ind in modelExcUps[3]:
    if (modelFlux.fluxes[ind] < 0):
      fluxLimit = fluxLimit + abs(modelFlux.fluxes[ind])
  # Making the "c" metabolite
  c = Metabolite('c', formula='', name='c')
  for ind in modelExcUps[3]:
    rxn = model.reactions[ind]
    # Skipping the "omit" reactions
    if (rxn.id in omit):
      continue
    # Adding "c" to the uptake reaction
    if (rxn.products == [] and rxn.upper_bound > 0):
      # Splitting bidirectional exchange reactions into uptake and export reactions
      rxnTMP = rxn.copy()
      rxnTMP.id = "TMP_" + rxnTMP.id
      rxnTMP.name = "TMP_" + rxnTMP.name
      rxnTMP.lower_bound = 0
      model.add_reaction(rxnTMP)
      rxn.reaction = (list(rxn.reactants)[0].id + " + c <=>")
      rxn.upper_bound = 0
    elif (rxn.products == [] and rxn.upper_bound == 0):
      rxn.reaction = (list(rxn.reactants)[0].id + " + c <=>")
  # Export - c
  reactionExpC = Reaction("EXP_c")
  reactionExpC.name = "EXP_c"
  reactionExpC.lower_bound = 0
  reactionExpC.upper_bound = fluxLimit
  reactionExpC.add_metabolites({c: -1.0})
  model.add_reaction(reactionExpC)
###########################################################################################

###########################################################################################
# A function for finding the high yielding media components
# Algorithm:
# 1. Implement 'c' to all the uptake reactions (e.g., x <=> to x + c <=>). Export reactions are excluded. Therefore, bidirectional
#    exchange reactions are split into uptake and export reactions. One can omit reaction(s) from adding 'c' by passing the reaction
#    IDs as an array via the "omit" parameter.
# 2. Implement a 'c' export reaction (c -->).
# 3. Ramp down 'c' from a high upper bound value (sum of all the uptake fluxes of the original model) to a value close to zero while
#  optimising for the BOF ("growthCutOff" parameter determines the minimum required growth rate). 
# 4. Count the negative flux carrying uptake reactions on each step ("numUptakeRxns").
# 5. Pick the uptake reactions of the lowest 'c' upper bound optimisation that also met the minimum growth rate requirement.
# These uptake reactions contain heigh-yield giving metabolites while optimising for the BOF.
def findHighYieldImport(model, omit, growthCutOff):
  j = 0
  maxUptake = 0
  maxUptakeFlux = 0
  step = 0
  numUptakeRxns = []
  addC(model, omit)
  modelExcUpsAddC = showExcRxns(model)
  maxUptake = len(modelExcUpsAddC[3])
  startValue = model.reactions.get_by_id("EXP_c").upper_bound
  # Identifying an appropriate step value
  step = 10**(len(str(int(startValue)))-1)
  while (1):
    if (startValue > 0):
      model.reactions.get_by_id("EXP_c").upper_bound = startValue
    else:
      model.reactions.get_by_id("EXP_c").upper_bound = 0
    bFluxTMP = model.optimize()
    if (bFluxTMP.objective_value > growthCutOff):
      # Quit if export restriction of the currency metabolite 'c' does not have any significant effect on the optima
      if (model.reactions.get_by_id("EXP_c").upper_bound == 0):
        print ("\nError: export restriction of the currency metabolite 'c' does not have any significant effect on the optima!\n")
        break
      # Counting the active uptake reactions
      for ind in modelExcUpsAddC[3]:
        if (bFluxTMP.fluxes[ind] < 0):
          j = j + 1
      if (j < maxUptake):
        maxUptake = j
        numUptakeRxns.append(j)
        maxUptakeFlux = bFluxTMP
      j = 0
      # Decreasing the "c" upper bound
      startValue = startValue - step
    else:
      startValue = startValue + step
      model.reactions.get_by_id("EXP_c").upper_bound = startValue
      if (startValue < 1): break
      step = step/10
      startValue = startValue - step    
  # Media composition
  mediaComposition = []
  mediaIndex = []
  print ("\nID and flux of the high-yield Minimal Medium uptake reactions:")
  for ind in modelExcUpsAddC[3]:
    if (maxUptakeFlux.fluxes[ind] < 0):
      mediaComposition.append(model.reactions[ind].id)
      print("{}\t{}".format(model.reactions[ind].id, maxUptakeFlux.fluxes[ind]))
  print ("\n")
  return mediaComposition, maxUptakeFlux, numUptakeRxns
###########################################################################################

###########################################################################################
# A function to find the necessary uptake reactions by turing one off at a time
def findCritical(model, mediaComposition, growthCutOff):
  growthRate = []
  impUptake = []
  for rxnId in mediaComposition:
    rxn = model.reactions.get_by_id(rxnId)
    rxn.lower_bound = 0
    growthRate.append(model.slim_optimize())
    if(growthRate[-1] < growthCutOff):
      impUptake.append(model.reactions.index(rxn))
    rxn.lower_bound = -1000
  return impUptake, growthRate
###########################################################################################

###########################################################################################
# A function to eliminate the non-critical uptake reactions that has the least effect on the growth rate at each iteration
def uptakeEliminator(model, upIndex, desirableGrowth):
  delInd = []
  i = len(upIndex)
  while i > 0:
    j = 0
    tDel = -1
    for ind in upIndex:
      if (ind not in delInd):
        rxn = model.reactions[ind]
        rxn.lower_bound = 0
        gRate = model.slim_optimize()
        if(gRate >= desirableGrowth):
          if (gRate > j):
            j = gRate
            tDel = ind
        rxn.lower_bound = -1000
    if (tDel == -1):
      break
    delInd.append(tDel)
    rxn = model.reactions[tDel]
    rxn.lower_bound = 0
    i = i - 1
  trimmedMedia = list(set(upIndex) - set(delInd))
  return trimmedMedia
###########################################################################################

###########################################################################################
# This funciton enables only the provided uptake reactions and checks if the model can be optimised at or above the growthCutOff.
# It then returns the flux carrying uptake reactions, flux values and the optimisation result. If 'display' is true, then prints the
# optimisation result, reactions and fluxes.
def findActiveUptake(model, upIndex, growthCutOff, display = 0):
  modelExcUps = []
  modelFlux = []
  fluxRxns = []
  uptakeFlux = []
  gRate = 0
  # Finding the uptake reactions
  modelExcUps = showExcRxns(model)
  # Restricting all the uptake reactions
  restrictBounds(model, modelExcUps[3], 0, 1000)
  # Enabling the provided uptake reactions and optimise
  restrictBounds(model, upIndex, -1000, 1000)
  modelFlux = model.optimize()
  gRate = modelFlux.objective_value
  if (gRate >= growthCutOff):
    if (display):
      print("\nOptimisation result: {}".format(gRate))
      print("\nUptake reactions and fluxes:")
    for ind in upIndex:
      if (modelFlux.fluxes[ind] < 0):
        fluxRxns.append(ind)
        uptakeFlux.append(modelFlux.fluxes[ind])
        if (display): print("{}\t{}".format(model.reactions[ind].id, modelFlux.fluxes[ind]))
    if (display): print("\n")
  return fluxRxns, uptakeFlux, gRate
###########################################################################################
