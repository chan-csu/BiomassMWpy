from __future__ import print_function
import cobra
#Add the path for importing the module "biomassmw"
#cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"BiomassMW","PythonCobrapy")))
#if cmd_subfolder not in sys.path:
#	sys.path.insert(0, cmd_subfolder)

from biomassmw import biomassmw

#The example model is the intermediate Bacteroides thetaiotaomicron iAH991 model after
#correcting the reactions (except proton imbalance, which can be identified
#by the procedure) but before renormalizing the biomass reaction.
model = cobra.io.load_matlab_model('example_model.mat')

for solver in ['gurobi']: #'glpk', 'cplex'
	model.solver = solver
	print("Solve using %s (object: %s)" %(solver, model.solver.__repr__()))
	balance = biomassmw(model)
	balance.compute_met_form(findCM='null',nameCM = 0)
	cm = balance.met_form_results.cm_info.cm_generic_dict
	for cmName, cmDict in  cm.items():
		print('Conserved moiety: %s' %cmName.formula)
		for i in cmDict:
			print('    %s\t%s' %(i.id, balance.met_form_results.formulae[i]))
	print(' ')
	model2 = balance.met_form_results.model
	bm = model2.metabolites[model2.metabolites.index('biomass[e]')]
	bm0 = balance.model.metabolites[balance.model.metabolites.index('biomass[e]')]
	print('Biomass formula: %s\nCharge: %.4f\nBiomass weight: %.4f g/mol' %(bm.formula, bm.charge, balance.met_form_results.formulae[bm0].mw))
	print('\nFind the range for biomass MW:')
	balance.compute_met_range('biomass[e]')
	print('[%.6f, %.6f]' %(balance.met_range_results.mw_range[0], balance.met_range_results.mw_range[1]))
