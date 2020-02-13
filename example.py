from __future__ import print_function
import cobra
from biomassmw import biomassmw

#The example model is the intermediate Bacteroides thetaiotaomicron iAH991 model after
#correcting the reactions (except proton imbalance, which can be identified
#by the procedure) but before renormalizing the biomass reaction.
model = cobra.io.load_matlab_model('example_model.mat')

for solver in ['gurobi']: #'glpk', 'cplex'
	model.solver = solver
	print("Solve using %s (object: %s)" %(solver, model.solver.__repr__()))
	# initialize an instance for balancing the model
	balance = biomassmw(model)
	# Functionality 1:
	# compute formulae for metabolites with unknown formulae
	# use null space to find the conserved quantities
	balance.compute_met_form(findCM='null', nameCM = 0)
	# print all the conserved quantities identified
	cm = balance.met_form_results.cm_info.cm_generic_dict
	for cmName, cmDict in  cm.items():
		print('Conserved moiety: %s' %cmName.formula)
		for i in cmDict:
			print('    %s\t%s' %(i.id, balance.met_form_results.formulae[i]))
	print(' ')
	# the model incorporated with newly computed formulae
	model2 = balance.met_form_results.model
	# biomass metabolite
	bm = model2.metabolites[model2.metabolites.index('biomass[e]')]
	bm0 = balance.model.metabolites[balance.model.metabolites.index('biomass[e]')]
	print('Biomass formula: %s\nCharge: %.4f\nBiomass weight: %.4f g/mol' %(bm.formula, bm.charge, balance.met_form_results.formulae[bm0].mw))
	# print unbalanced reactions
	print('\nUnbalanced reactions:\n')
	nUnbalRxn = 0
	for rxn, bal in balance.met_form_results.rxn_balance.items():
		if bal != {} and len(rxn.metabolites) > 1:
			nUnbalRxn += 1
			print("%-2d: %20s | %s" %(nUnbalRxn, rxn.id, ', '.join([e + ": %f" %(b) for e,b in bal.items() ])))
	# most of the unbalanced reactions have just proton imbalance.

	# Functionality 2:
	# compute the range of MW for the biomass metabolite assuming under the condition of minimum inconsistency in elemental balance
	print('\nFind the range for biomass MW:')
	balance.compute_met_range('biomass[e]')
	print('[%.6f, %.6f]' %(balance.met_range_results.mw_range[0], balance.met_range_results.mw_range[1]))
