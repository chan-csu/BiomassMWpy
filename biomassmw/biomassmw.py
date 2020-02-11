from __future__ import print_function
from itertools import chain
import re
import pdb
from cobra.core import Model, Reaction, Metabolite
from cobra.util.solver import set_objective
from .objects.generic_formula import generic_formula as Formula
from .functions import solution_infeasibility, active_met_rxn, formula_dict2str, num2alpha, extreme_rays_from_null_basis
from .objects.result_structure import data_object, preprocess_data, min_incon_parsi_info, conserved_moiety_info, result_structure
import datetime
try:
	import cdd
	cddImported = True
except ImportError:
	print('pycddlib is not installed. Unable to use it for extreme ray calculations. Calculate the null space biasis only')
	cddImported = False

element_re = re.compile("([A-Z][a-z_]*)(\-?[0-9.]+[0-9.]?|(?=[A-Z])?)")

class biomassmw(data_object):
	'''biomassmw is a class for performing minimum inconsistency under parsimony (MIP) to determine the chemical formulae or the range for the molecular weight of a metabolite
	Must be initialized with a cobra model: biomassmw(cobra_model)
	'''

	#this class property determines whether conserved moiety calculation can be performed
	__cddImported = cddImported

	def __init__(self, cobra_model):
		#must be associated with a normal cobra model
		if isinstance(cobra_model, Model):
			self.model = cobra_model
		else:
			raise ValueError('biomassmw must be initiated with a cobrapy model')
		#default maximum value for variables
		self.infinity = float('inf')
		#default minimum value for charge variables
		self.negative_infinity = - 10000 # -float('inf') seems to cause error for some solvers
		self.__solver = cobra_model.solver

	def compute_met_form(self, met_known=None, rxns=None, met_fill=['HCharge'], findCM='null', deadend=True, nameCM=0, **kwargs):
		'''Compute the chemical formulas of the unknown metabolites given the formulae for a set of known metabolites using a set of reactions.
		met_known: metabolites with known formulae, in a list or Dictlist of metabolite objects or a list of metabolite IDs. Defaulted all metabolites with formulae.
		rxns: reactions used for balancing, in a list or Dictlist of metabolite objects or a list of metabolite IDs. Defaulted all active non-exchange reactions.
		met_fill: metabolites for automatically filling inconsistency, in a list of chemical formulae. Default 'HCharge' (proton)
		findCM: Method used to calculate extreme rays for conserved moieties. 'null' (default, rational null basis from reduced row echelon form) or 'cdd' (using the package pycddlib). Skip the calculation for all other values.
		deadend: Include deadend metabolites or not when calculating extreme rays (default True, include)
		nameCM: 0 to use defaulted names 'Conserve_xx' for all moieties.
		        1 to use defaulted names only for conserved moieties involved in deadend metabolites. Interactively name others.
		        2 to name all conserved moieties identified.
		**kwargs: name-value parameters for optimizing cobrapy models, e.g. tolerance_feasibility=1e-9
		Return a ResultStructure object with the following attributes:
			formulae: the final formulae for unknonw metabolites
			S_fill: a dictionary of stoichiometric coefficients for filling metabolites
			rxn_balance: {rxn:{e: balance}} dictionary for the elemental balance for each element e of each reaction rxn
			model: a copy of the input cobra model with updated chemical formulae
			mip_info: min_incon_parsi_info object summarizing the results of solving MIP.
			cm_info: conserved_moiety_info object containing the conserved moiety information.
		Also avaialble as self.met_form_results
		'''
		if not hasattr(self, 'pre'):
			#Pre-porocessing
			print('Preprocessing ... %s' %datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
			self.preprocessing(met_known, rxns, met_fill, **kwargs)
		#Solve 'Minimum Inconsistency under Parsimony' to find a set of minimal formulae
		print('Solving Minimum Inconsistency under Parsimony ... %s' %datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
		mip_info = self.min_incon_parsi(**kwargs)
		#Find conserved moieties by calculating extreme rays
		if not self.__cddImported and findCM =='cdd':
			print('pycddlib is not installed. Use rational null space basis.')
			findCM = 'null'
		cm_info = None
		if findCM in ['null', 'cdd'] and mip_info.final != 'infeasible':
			print('Finding conserved moieties ... %s' %datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
			cm_info = self.conserved_moieties(mip_info, findCM, deadend, nameCM)
		#get the updated model and info on solutions
		met_form_results = self.get_results_from_mip_and_cm(mip_info, cm_info)
		self.met_form_results = met_form_results
		return met_form_results

	def compute_met_range(self, metInterest, met_known=None, rxns=None, **kwargs):
		'''Compute the minimum and maximum possible MW of the target metabolite
		metInterest: metabolite of interest, in a metabolite object or ID.
		met_known: metabolites with known formulae, in a list or Dictlist of metabolite objects or a list of metabolite IDs. Defaulted all metabolites with formulae.
		rxns: reactions used for balancing, in a list or Dictlist of metabolite objects or a list of metabolite IDs. Defaulted all active non-exchange reactions.
		**kwargs: name-value parameters for optimizing cobrapy models, e.g. tolerance_feasibility=1e-9
		Return min_incon_parsi_info object summarizing the results of solving MIP.
			formulae: (formula for min MW, formula max MW)
			mw_range: (min MW, max MW)
			rhs: {e: {i: RHS[i] for all metabolite i}} the RHS value in the MIP problem for each element solved. (compute_met_range only)
			infeas: infeasibility of each solve
			bound: bound used for total inconsistency or the relaxation value eps0 for each solve
			obj: objective function value for each solve
			solution: solution values for each type of variables (m, xp, xn)
			met_model: the MIP problem solved for each element as a cobra model. Same model but different rhs for different elements.
			final: the final status of the solution
			sol_stat: solution status for each element/connected set of elements
		Also avaialble as self.met_range_results
		'''
		if not hasattr(self, 'pre'):
			#Pre-porocessing
			print('Preprocessing ... %s' %datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
			self.preprocessing(met_known, rxns, **kwargs)
		#Solve 'Minimum Inconsistency under Parsimony' to find a set of minimal formulae
		mip_info = self.min_incon_parsi_mwRange(metInterest, **kwargs)
		self.met_range_results = mip_info
		return mip_info

	def preprocessing(self, met_known=None, rxns=None, met_fill = ['HCharge'], **kwargs):
		'''Preprocessing for running min_incon_parsi or min_incon_parsi_mwRange. Called by compute_met_form or compute_met_range
		'''
		model = self.model
		#metabolites with known formulae
		if met_known is None:
			#default all mets with non-generic chemical formulae
			metK = {i: Formula(i) for i in model.metabolites}
			metK = {i: i2 for i, i2 in metK.items() if not i2.unknown} #and not i2.generic}
		else:
			if (type(met_known) is list or isinstance(met_known,type(model.metabolites)))\
			and all([isinstance(i, type(model.metabolites[0])) for i in met_known]):
				#Dictlist or list of met objects
				metK = {i: Formula(i) for i in met_known if i in model.metabolites}
			elif all([isinstance(i, str) for i in met_known]):
				#list of met IDs
				metK = {model.metabolites.get_by_id(i): Formula(model.metabolites[model.metabolites.index(i)]) for i in met_known if i in model.metabolites}
			else:
				raise	ValueError('met_known must be a list of met IDs or met objects.')
			if len(metK) < len(met_known):
				raise	ValueError('Some mets in met_known are not in the model')
			if not all([bool(i.formula) for i in metK]):
				print('Some mets in met_known do not have formula. Ignore them.')
				metK = {i:i2 for i,i2 in metK.items() if bool(i2.formula)}
		#unknown mets
		metU = [i for i in model.metabolites if i not in metK]
		#reactions to be balanced for minimum inconsistency
		if rxns is None:
			#all non-exchange reactions
			#rxns = {model.reactions.index(j):j for j in model.reactions if bool(j.reactants) and bool(j.products)}
			rxnK = [j for j in model.reactions if bool(j.reactants) and bool(j.products) and (j.upper_bound != 0 or j.lower_bound != 0)]
		else:
			if (type(rxns) is list or isinstance(rxns,type(model.reactions))) \
			and all([isinstance(j, Reaction) for j in rxns]):
				#Dictlist of rxn objects
				rxnK = [j for j in rxns if j in model.reactions]
			elif all([isinstance(j, str) for j in rxns]):
				#list of rxn IDs
				rxnK = [model.reactions[model.reactions.index(j)] for j in rxns if j in model.reactions]
			else:
				raise ValueError('rxns must be a list of rxn IDs or rxn objects.')
			if len(rxnK) < len(rxns):
				raise ValueError('Some rxns are not in the model')

		#metabolites for filling inconsistency
		metF = {met_fill[i]: Formula(met_fill[i]) for i in range(len(met_fill))} if met_fill is not None else None

		#find elements that are connected because of met_fill. They need to be optimized in the same problem.
		ele_connect, met_fill_connect = [], {}
		eleSet = set([e for i in chain(metK.values(), metF.values()) for e in i.elements])
		ele = list(eleSet)
		while bool(eleSet):
			eleCur = set([eleSet.pop()])
			while True:
				met_fillCur = [i for i in metF if bool(eleCur.intersection(metF[i].elements))]
				eleNext = set([e for i in met_fillCur for e in metF[i].elements]).union(eleCur) #in case of met_fillCur = []
				if eleCur == eleNext:
					break
				eleCur = eleNext.copy()
			ele_connect.append(tuple(eleCur))
			met_fill_connect[tuple(eleCur)] = met_fillCur
			eleSet.difference_update(eleCur)

		#feasibility tolerance
		feasTol = kwargs['tolerance_feasibility'] if 'tolerance_feasibility' in kwargs else 1e-8
		#digits rounded to for avoiding numerical issues
		digitRounded = 10

		pre = preprocess_data()
		pre.met_known, pre.met_unknown, pre.met_fill, pre.rxn_known, pre.ele, pre.ele_connect, pre.met_fill_connect \
		,pre.feasTol, pre.digitRounded = metK, metU, metF, rxnK, ele, ele_connect, met_fill_connect, feasTol, digitRounded
		self.pre = pre
		return

	def min_incon_parsi(self, **kwargs):
		'''The main step called by compute_met_form to solve the MIP problem.
		Return a min_incon_parsi_info object summarizing the results of solving MIP.
			formulae: result formulae
			infeas: infeasibility of each solve
			bound: bound used for total inconsistency or the relaxation value eps0 for each solve
			obj: objective function value for each solve
			solution: solution values for each type of variables (m, xp, xn, Ap, An)
			met_model: the MIP problems solved for each connected set of elements as individual cobra models
			final: the final status of the solution
			sol_stat: solution status for each connected set of elements
			sol_constrain: 'minIncon' or 'minFill' indicating the solution used to constrain the minimal formulae problem.
		'''
		inf, neg_inf = self.infinity, self.negative_infinity
		pre = self.pre
		metK, metU, metF, rxnK, ele, ele_connect, met_fill_connect, feasTol, digitRounded \
		= pre.met_known, pre.met_unknown, pre.met_fill, pre.rxn_known, \
		pre.ele, pre.ele_connect, pre.met_fill_connect, pre.feasTol, pre.digitRounded
		model = self.model
		kwargs['objective_sense'] = 'minimize' #always minimize
		#data to be stored
		infeas, bound, solution, met_model, sol_constrain, sol_stat, obj = ({} for i in range(7))
		for k in ['minIncon', 'minFill', 'minForm']: #pre-assignment
			solution[k] = {'m': {}, 'xp': {}, 'xn': {}, 'Ap': {}, 'An': {}}

		#Optimize for each connected componenet in ele_connect
		ct = 0
		for eCC in ele_connect:
			ct += 1
			print("Optimizing for %d / %d set of connected elements ... %s" %(ct, len(ele_connect), datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
			metModelJ = Model(', '.join(eCC))
			metModelJ.solver = self.__solver
			infeasJ, boundJ, objJ = {}, {}, {}
			constraint = {e: {j: Metabolite(j.id + ',' + e) for j in rxnK} for e in eCC}
			m = {e: {i: Reaction('m_' + i.id + ',' + e) for i in metU} for e in eCC}
			xp = {e: {j: Reaction('xp_' + j.id + ',' + e) for j in rxnK} for e in eCC}
			xn = {e: {j: Reaction('xn_' + j.id + ',' + e) for j in rxnK} for e in eCC}
			#add adjustment variable if filling mets associated with the current connected elements
			Ap = {i: {j: Reaction('Ap_' + metF[i].formula + ',' + j.id) for j in rxnK} for i in met_fill_connect[eCC]}
			An = {i: {j: Reaction('An_' + metF[i].formula + ',' + j.id) for j in rxnK} for i in met_fill_connect[eCC]}
			for e in eCC:
				for j in rxnK:
					#RHS for each constraint: -sum(S_ij * m^known_ie) (obsolete)
					constraint[e][j]._bound = -sum([S_ij * metK[i].elements[e] for i, S_ij in j._metabolites.items() if i in metK and e in metK[i].elements])
					constraint[e][j]._constraint_sense = 'E'
					#x_pos - x_neg for each constraint
					xp[e][j].add_metabolites({constraint[e][j]: 1})
					xn[e][j].add_metabolites({constraint[e][j]: -1})
					xp[e][j].lower_bound, xn[e][j].lower_bound, xp[e][j].upper_bound, xn[e][j].upper_bound = 0, 0, inf, inf
					#m * A_pos - m * A_neg for each constraint
					for i in met_fill_connect[eCC]:
						if e in metF[i].elements:
							Ap[i][j].add_metabolites({constraint[e][j]:  metF[i].elements[e]})
							An[i][j].add_metabolites({constraint[e][j]: -metF[i].elements[e]})
				for i in metU:
					# S_ij x m_ie for each i
					m[e][i].add_metabolites({constraint[e][j]: j._metabolites[i] for j in list(i._reaction) if j in rxnK})
					m[e][i].upper_bound = inf
					m[e][i].lower_bound = neg_inf if e == 'Charge' else 0
			for i in met_fill_connect[eCC]:
				for j in rxnK:
					Ap[i][j].lower_bound, An[i][j].lower_bound, Ap[i][j].upper_bound, An[i][j].upper_bound = 0, 0, inf, inf

			#add reactions into the model
			metModelJ.add_reactions([m[e][i] for e in eCC for i in metU])
			metModelJ.add_reactions([xp[e][j] for e in eCC for j in rxnK])
			metModelJ.add_reactions([xn[e][j] for e in eCC for j in rxnK])
			if bool(met_fill_connect[eCC]):
				metModelJ.add_reactions([Ap[i][j] for i in met_fill_connect[eCC] for j in rxnK])
				metModelJ.add_reactions([An[i][j] for i in met_fill_connect[eCC] for j in rxnK])
            # the objective can be set as a dictionary of rxns {rxn_id: 1,...}
			objective_dict = {xp[e][j]: 1 for e in eCC for j in rxnK}
			objective_dict.update({xn[e][j]: 1 for e in eCC for j in rxnK})
			set_objective(metModelJ, objective_dict)
			#set RHS: -sum(S_ij * m^known_ie)
			for e in eCC:
				for j in rxnK:
					metModelJ.constraints[constraint[e][j].id].ub = inf #ub must be set to be > lb to avoid error
					metModelJ.constraints[constraint[e][j].id].lb, metModelJ.constraints[constraint[e][j].id].ub = constraint[e][j]._bound, constraint[e][j]._bound

			#Solve for minimum inconsistency
			sol = metModelJ.optimize(**kwargs)
			solStatJ = 'minIncon'
			infeasJ[solStatJ] = solution_infeasibility(metModelJ, sol)
			if sol.fluxes is None:
				boundJ[solStatJ] = {e: float('nan') for e in eCC}
			else:
				boundJ[solStatJ] = {e: sum([sol.fluxes[j.id] for j in chain(xp[e].values(), xn[e].values())]) for e in eCC}
			if not infeasJ[solStatJ] <= feasTol:
				#infeasible (should not happen)
				infeasJ['minFill'], infeasJ['minForm'] = inf, inf
				solConstrainJ = 'infeasible'
				solStatJ = 'infeasible'
				objJ['minIncon'], objJ['minFill'], objJ['minForm'] = (float('nan') for i in range(3))
			else:
				#Feasible. Store the solution
				solution[solStatJ]['m'].update({e: {i: sol.fluxes[m[e][i].id] for i in metU} for e in eCC})
				solution[solStatJ]['xp'].update({e: {j: sol.fluxes[xp[e][j].id] for j in rxnK} for e in eCC})
				solution[solStatJ]['xn'].update({e: {j: sol.fluxes[xn[e][j].id] for j in rxnK} for e in eCC})
				solution[solStatJ]['Ap'].update({i: {j: sol.fluxes[Ap[i][j].id] for j in rxnK} for i in met_fill_connect[eCC]})
				solution[solStatJ]['An'].update({i: {j: sol.fluxes[An[i][j].id] for j in rxnK} for i in met_fill_connect[eCC]})
				objJ[solStatJ] = sol.objective_value  # .f deprecated since 2018 Oct
				#solution used to constrain the minimal formula problem
				solConstrainJ = 'minIncon'
				#minimize total adjustment if filling metabolites exist
				if bool(met_fill_connect[eCC]):
					#Add constraint to fix the total inconsistency for each element
					constraint_minIncon = {e: Metabolite('minIncon_'+e) for e in eCC}
					for e in eCC:
						#rounding to avoid infeasibility due to numerical issues
						constraint_minIncon[e]._bound = round(boundJ[solStatJ][e], digitRounded)
						constraint_minIncon[e]._constraint_sense = 'L'
						#sum(xp) + sum(xn) <= total inconsistency
						for j in rxnK:
							xp[e][j].add_metabolites({constraint_minIncon[e]: 1})
							xn[e][j].add_metabolites({constraint_minIncon[e]: 1})

						metModelJ.constraints[constraint_minIncon[e].id].lb = neg_inf
						metModelJ.constraints[constraint_minIncon[e].id].ub = round(boundJ[solStatJ][e], digitRounded)

					#reset the objective function
					objective_dict = {Ap[i][j]: 1 for i in met_fill_connect[eCC] for j in rxnK}
					objective_dict.update({An[i][j]: 1 for i in met_fill_connect[eCC] for j in rxnK})
					set_objective(metModelJ, objective_dict)
					solStatJ = 'minFill'
					eps0 = 1e-6
					while True:
						sol = metModelJ.optimize(**kwargs)
						infeasJ[solStatJ] = solution_infeasibility(metModelJ, sol)
						if infeasJ[solStatJ] <= feasTol or eps0 > 1e-4 + 1e-8:
							break
						eps0 *= 10
						for e in eCC:
							#rounding to avoid infeasibility due to numerical issues
							metModelJ.constraints[constraint_minIncon[e].id].ub = round(boundJ['minIncon'][e] * (1 + eps0), digitRounded)

					boundJ[solStatJ] = eps0
					if infeasJ[solStatJ] <= feasTol:
						#Feasible. Use this as the solution for constraining the minimal formula problem
						solConstrainJ = 'minFill'
						#Store the solution
						solution[solStatJ]['m'].update({e: {i: sol.fluxes[m[e][i].id] for i in metU} for e in eCC})
						solution[solStatJ]['xp'].update({e: {j: sol.fluxes[xp[e][j].id] for j in rxnK} for e in eCC})
						solution[solStatJ]['xn'].update({e: {j: sol.fluxes[xn[e][j].id] for j in rxnK} for e in eCC})
						solution[solStatJ]['Ap'].update({i: {j: sol.fluxes[Ap[i][j].id] for j in rxnK} for i in met_fill_connect[eCC]})
						solution[solStatJ]['An'].update({i: {j: sol.fluxes[An[i][j].id] for j in rxnK} for i in met_fill_connect[eCC]})
						objJ[solStatJ] = sol.objective_value
					else:
						#infeasible, should not happen
						objJ[solStatJ] = float('nan')

					#prepare to compute minimal formulae
					eps0 = 1e-10
					for j in rxnK:
						for i in met_fill_connect[eCC]:
							#Fix the adjustment
							Ap[i][j].lower_bound = round(solution[solConstrainJ]['Ap'][i][j] * (1 - eps0), digitRounded)
							Ap[i][j].upper_bound = round(solution[solConstrainJ]['Ap'][i][j] * (1 + eps0), digitRounded)
							An[i][j].lower_bound = round(solution[solConstrainJ]['An'][i][j] * (1 - eps0), digitRounded)
							An[i][j].upper_bound = round(solution[solConstrainJ]['An'][i][j] * (1 + eps0), digitRounded)
						for e in eCC:
							#remove all coefficients on the constraints for total inconsistency
							tmp = xp[e][j]._metabolites.pop(constraint_minIncon[e])
							tmp = xn[e][j]._metabolites.pop(constraint_minIncon[e])
					for e in eCC:
						#then remove the constraints
						metModelJ.metabolites.remove(constraint_minIncon[e])
				else:
					#no filling metabolites
					infeasJ['minFill'], boundJ['minFill'] = 0, 0

				#compute minimal formulae
				eps0 = 1e-10
				for e in eCC:
					for j in rxnK:
						#fix the inconsistency
						xp[e][j].lower_bound = round(solution[solConstrainJ]['xp'][e][j] * (1 - eps0), digitRounded)
						xp[e][j].upper_bound = round(solution[solConstrainJ]['xp'][e][j] * (1 + eps0), digitRounded)
						xn[e][j].lower_bound = round(solution[solConstrainJ]['xn'][e][j] * (1 - eps0), digitRounded)
						xn[e][j].upper_bound = round(solution[solConstrainJ]['xn'][e][j] * (1 + eps0), digitRounded)

					if e == 'Charge':
						#add variables for the positive and negative part of charges
						chargePos = {i: Reaction('chargePos_' + i.id) for i in metU}
						chargeNeg = {i: Reaction('chargeNeg_' + i.id) for i in metU}
						constraint_chargeDecomp = {i: Metabolite('chargeDecomp_' + i.id) for i in metU}
						for i in metU:
							#constraint_chargeDecomp[i]: m_charge,i - chargePos_i + chargeNeg_i = 0
							m[e][i].add_metabolites({constraint_chargeDecomp[i]: 1})
							chargePos[i].add_metabolites({constraint_chargeDecomp[i]: -1})
							chargeNeg[i].add_metabolites({constraint_chargeDecomp[i]: 1})
							chargePos[i].lower_bound, chargePos[i].upper_bound = 0, inf
							chargeNeg[i].lower_bound, chargeNeg[i].upper_bound = 0, inf
							constraint_chargeDecomp[i]._bound, constraint_chargeDecomp[i]._constraint_sense = 0, 'E'
						metModelJ.add_metabolites(constraint_chargeDecomp.values())
						metModelJ.add_reactions(list(chargePos.values()) + list(chargeNeg.values()))

				#reset the objective function
				objective_dict = {m[e][i]: 1 for e in eCC if e != 'Charge' for i in metU}
				if 'Charge' in eCC:
					objective_dict.update({chargePos[i]: 1 for i in metU})
					objective_dict.update({chargeNeg[i]: 1 for i in metU})
				set_objective(metModelJ, objective_dict)

				#Solve for minimum formulae
				solStatPrev, solStatJ =solStatJ, 'minForm'
				while True:
					sol = metModelJ.optimize(**kwargs)
					infeasJ[solStatJ] = solution_infeasibility(metModelJ, sol)
					if infeasJ[solStatJ] <= feasTol or eps0 > 1e-5 + 1e-8:
						break
					eps0 *= 10
					#relax bounds, rounding to avoid infeasibility due to numerical issues
					for j in rxnK:
						for i in met_fill_connect[eCC]:
							Ap[i][j].lower_bound = round(solution[solConstrainJ]['Ap'][i][j] * (1 - eps0), digitRounded)
							Ap[i][j].upper_bound = round(solution[solConstrainJ]['Ap'][i][j] * (1 + eps0), digitRounded)
							An[i][j].lower_bound = round(solution[solConstrainJ]['An'][i][j] * (1 - eps0), digitRounded)
							An[i][j].upper_bound = round(solution[solConstrainJ]['An'][i][j] * (1 + eps0), digitRounded)
						for e in eCC:
							xp[e][j].lower_bound = round(solution[solConstrainJ]['xp'][e][j] * (1 - eps0), digitRounded)
							xp[e][j].upper_bound = round(solution[solConstrainJ]['xp'][e][j] * (1 + eps0), digitRounded)
							xn[e][j].lower_bound = round(solution[solConstrainJ]['xn'][e][j] * (1 - eps0), digitRounded)
							xn[e][j].upper_bound = round(solution[solConstrainJ]['xn'][e][j] * (1 + eps0), digitRounded)
				boundJ[solStatJ] = eps0
				if infeasJ[solStatJ] <= feasTol:
					#Feasible. Store the solution
					solution[solStatJ]['m'].update({e: {i: sol.fluxes[m[e][i].id] for i in metU} for e in eCC})
					solution[solStatJ]['xp'].update({e: {j: sol.fluxes[xp[e][j].id] for j in rxnK} for e in eCC})
					solution[solStatJ]['xn'].update({e: {j: sol.fluxes[xn[e][j].id] for j in rxnK} for e in eCC})
					solution[solStatJ]['Ap'].update({i: {j: sol.fluxes[Ap[i][j].id] for j in rxnK} for i in met_fill_connect[eCC]})
					solution[solStatJ]['An'].update({i: {j: sol.fluxes[An[i][j].id] for j in rxnK} for i in met_fill_connect[eCC]})
					objJ[solStatJ] = sol.objective_value
				else:
					#infeasible, should not happen
					solStatJ = solStatPrev
					objJ[solStatJ] = float('nan')

			#store data
			infeas[eCC] = infeasJ
			bound[eCC] = boundJ
			obj[eCC] = objJ
			met_model[eCC] = metModelJ
			sol_constrain[eCC] = solConstrainJ
			sol_stat[eCC] = solStatJ

		#summarize the final solution state
		if any([k == 'infeasible' for k in sol_stat.values()]):
			print('Failure: no feasible solution can be found.')
			solFinal = 'infeasible'
		else:
			solMixed = True
			for stat in ['minForm', 'minFill', 'minIncon']:
				if all([k == stat for k in sol_stat.values()]):
					solFinal, solMixed = stat, False
					break
			if solMixed:
				solFinal = 'mixed'

		#Get the resultant set of formulae. For each set of elements in ele_connect, choose the latest solution (minForm > minFill > minIncon), recorded in sol_stat.
		if solFinal != 'infeasible':
			formulae = {i: Formula('Mass0') for i in metU}
			for i in metU:
				for eCC in ele_connect:
					formulae[i].update_elements({e: round(solution[sol_stat[eCC]]['m'][e][i], digitRounded) for e in eCC})
		else:
			formulae = {i: Formula() for i in metU}

		formulae.update(metK)
		mip_info = min_incon_parsi_info()
		mip_info.formulae, mip_info.infeas, mip_info.bound, mip_info.obj, mip_info.solution, \
		mip_info.met_model, mip_info.final, mip_info.sol_constrain, mip_info.sol_stat \
		= formulae, infeas, bound, obj, solution, met_model, solFinal, sol_constrain, sol_stat

		return mip_info

	def conserved_moieties(self, mip_info, findCM='null', deadend=True, nameCM=0):
		'''Find conserved moieties by computing extreme rays. Called by compute_met_form
		Return conserved_moiety_info object containing the conserved moiety information:
		cm: list of conserved moieties, each being a metabolite-coefficient dictionary
		cm_generic: list of entries in cm that are generic (no known metabolites involved)
		cm_generic_dict: {fomrula: cm_generic} dictionary. formula is the generic_formula object with the defaulted or inputted name as the formula for the conserved moiety.
		'''
		model = self.model
		if (not deadend) or nameCM == 1:
			activeMets, activeRxns = active_met_rxn(model)
		else:
			activeMets, activeRxns = model.metabolites, model.reactions
		if not deadend:
			cmMets, cmRxns = activeMets, activeRxns
		else:
			cmMets, cmRxns = model.metabolites, model.reactions

		#(matrix format: [[b_1, a_11, a_12, ..., a_1N], ..., [b_M, a_M1, a_M2, ..., a_MN]] for Ax + b >= 0
		#where A = [a_ij], b = [b_1, ..., b_M])

		if findCM == 'null':
			#transpose of S
			S = [[j._metabolites[i] if i in j._metabolites else 0 for i in cmMets] for j in cmRxns]
			#This method calculates a rational basis for transpose(S) from the reduced row echelon form, usually find a subset of extreme rays, quite probably the whole set.
			N = extreme_rays_from_null_basis(S)
			cm = [{cmMets[i]: N[i,k] for i in range(len(cmMets)) if N[i,k] != 0} for k in range(N.shape[1])]
		elif findCM == 'cdd':
			#transpose(S) >= 0
			S = [[0] + [j._metabolites[i] if i in j._metabolites else 0 for i in cmMets] for j in cmRxns]
			# #transpose(S) <= 0
			S += [[-i for i in j] for j in S]
			#all entries >= 0
			S += [[0] + [1 if i == j else 0 for i in range(len(cmMets))] for j in range(len(cmMets))]
			print('Matrix size for cdd extreme ray calculation: %d x %d' %(len(S), len(S[0])))
			# The cdd library seems unable to cope with genome-scale models. The best is to call EFMtool. To be implemented.
			mat = cdd.Matrix(S, number_type='float')
			mat.rep_type = cdd.RepType.INEQUALITY
			poly = cdd.Polyhedron(mat)
			ext = poly.get_generators()
			cm = [{cmMets[i]: ext.__getitem__(k)[i+1] for i in range(len(cmMets)) if ext.__getitem__(k)[i+1] != 0} for k in range(ext.row_size)]

		#generic conserved moieties involing no known metabolites
		cm_generic = [c for c in cm if all([i in self.pre.met_unknown for i in c])]
		cm_generic_dict = {}
		NcmDefault = 0
		for c in cm_generic:
			#Use defaulted names for dead end metabolites if nameCM = 1, or always use defaulted names if nameCM = 0
			if nameCM == 0 or (any([not i in activeMets for i in c.keys()]) and nameCM == 1):
				#defaulted names
				NcmDefault += 1
				cmNameCur = 'Conserve_' + num2alpha(NcmDefault)
			else:
				print('\n\n')
				for i in c.keys():
					toPrint = self.pre.met_known[i].formula if i in self.pre.met_known else mip_info.formulae[i].formula
					if toPrint == 'Mass0':
						toPrint = ''
					toPrint += formula_dict2str({"Conserve": c[i]})
					print('%s\t%s\t%s' %(i.id, i.name, toPrint))
				while True:
					cmNameCur = raw_input("\nEnter the formula for the conserved moiety: " \
					+ "(e.g. C2HRab_cd0.5Charge-1 -> {C:2, H:1, Rab_cd: 0.5, Charge: -1}, " \
					+ "hit return to use default name 'Conserve_xxx')\n")
					#check if the input is empty or a correct formula
					if cmNameCur == "" or ''.join([''.join(k) for k in element_re.findall(cmNameCur)]) == cmNameCur:
						break
					print('Incorrect format of the input formula!\n')
				if cmNameCur == '':
					#empty string means using the default name
					NcmDefault += 1
					cmNameCur = 'Conserve_' + num2alpha(NcmDefault)

			cm_generic_dict[Formula(cmNameCur)] = c

		cm_info = conserved_moiety_info()
		cm_info.cm, cm_info.cm_generic, cm_info.cm_generic_dict = cm, cm_generic, cm_generic_dict
		return cm_info

	def get_results_from_mip_and_cm(self, mip_info, cm_info=None):
		'''Summarize the final results. Called by compute_met_form
		'''
		met_form_results = result_structure()

		if mip_info.final =='infeasible':
			met_form_results.mip_info = mip_info
			met_form_results.model, met_form_results.formulae, met_form_results.S_fill, met_form_results.rxn_balance, \
			met_form_results.mip_info, met_form_results.cm_info = (None for k in range(6))
			self.met_form_results = met_form_results
			return met_form_results

		#Get S_fill. For each set of elements in ele_connect, choose the latest solution (minForm > minFill > minIncon), recorded in mip_info.sol_stat.
		S_fill = {}
		for eCC in self.pre.ele_connect:
			for i in self.pre.met_fill_connect[eCC]:
				d = {j: round(mip_info.solution[mip_info.sol_stat[eCC]]['Ap'][i][j] - mip_info.solution[mip_info.sol_stat[eCC]]['An'][i][j], self.pre.digitRounded) for j in self.pre.rxn_known
					if abs(mip_info.solution[mip_info.sol_stat[eCC]]['Ap'][i][j] - mip_info.solution[mip_info.sol_stat[eCC]]['An'][i][j]) > 1e-6}
				if bool(d):
					S_fill[i] = d

		formulae = mip_info.formulae.copy()
		if cm_info != None and bool(cm_info.cm_generic_dict):
			#include the unknown conserved moieties into the chemical formulae
			for cmNameCur, c in cm_info.cm_generic_dict.items():
				for i in c:
					formulae[i].update_elements({e: c[i]*m_ie for e, m_ie in cmNameCur.elements.items()}, True)

		#elemental balance of each reaction
		rxn_balance = {j: {e: sum([j._metabolites[i] * formulae[i].elements[e] for i in j._metabolites if e in formulae[i].elements])
			for e in self.pre.ele} for j in self.model.reactions}
		rxn_balance = {j: {e: rxn_balance[j][e] for e in rxn_balance[j] if rxn_balance[j][e] != 0} for j in rxn_balance}
		#updated model
		model = self.model.copy()
		for i in model.metabolites:
			#the original metabolite object
			i0 = self.model.metabolites[self.model.metabolites.index(i.id)]
			i.formula, i.charge = formulae[i0].formula_wo_charge, formulae[i0].charge

		met_form_results.model, met_form_results.formulae, met_form_results.S_fill, met_form_results.rxn_balance, \
		met_form_results.mip_info, met_form_results.cm_info = model, formulae, S_fill, rxn_balance, mip_info, cm_info
		self.met_form_results = met_form_results
		return met_form_results

	def min_incon_parsi_mwRange(self, metInterest, **kwargs):
		'''The main step called by compute_met_range to find the range for the molecular weight
		Return min_incon_parsi_info object summarizing the results of solving MIP.
			formulae: (formula for min MW, formula max MW)
			mw_range: (min MW, max MW)
			rhs: {e: {i: RHS[i] for all metabolite i}} the RHS value in the MIP problem for each element solved. (compute_met_range only)
			infeas: infeasibility of each solve
			bound: bound used for total inconsistency or the relaxation value eps0 for each solve
			obj: objective function value for each solve
			solution: solution values for each type of variables (m, xp, xn)
			met_model: the MIP problem solved for each element as a cobra model. Same model but different rhs for different elements.
			final: the final status of the solution
			sol_stat: solution status for each element/connected set of elements
		'''
		inf, neg_inf = self.infinity, self.negative_infinity
		pre = self.pre
		metK, metU, rxnK, ele, feasTol, digitRounded \
		= pre.met_known, pre.met_unknown, pre.rxn_known, pre.ele, pre.feasTol, pre.digitRounded
		model = self.model
		#handle metInterest
		if isinstance(metInterest,type(model.metabolites[0])):
			metI = metInterest
		elif isinstance(metInterest, str):
			metI = model.metabolites.get_by_id(metInterest)
		if metI in metK:
			print("%s in the input is already known." %metI.id)
			mip_info = min_incon_parsi_info()
			mip_info.mw_range = (metK[metI].mw, metK[metI].mw)
			return mip_info

		print('Find the range for the molecular weight of %s ... %s' %(metI.id, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

		infeas, bound, solution, sol_stat, obj, rhs = ({} for i in range(6))
		for k in ['minIncon', 'minMw', 'maxMw']: #pre-assignment
			solution[k] = {'m': {}, 'xp': {}, 'xn': {}}
		ct = 0
		#The optimization problem for each element is the same except the RHS
		met_model = Model('min/max Mw of ' + metI.id)
		met_model.solver = self.__solver
		constraint = {j: Metabolite(j.id) for j in rxnK}
		m = {i: Reaction('m_' + i.id) for i in metU}
		xp = {j: Reaction('xp_' + j.id) for j in rxnK}
		xn = {j: Reaction('xn_' + j.id) for j in rxnK}
		for j in rxnK:
			xp[j].add_metabolites({constraint[j]: 1})
			xn[j].add_metabolites({constraint[j]: -1})
			xp[j].lower_bound, xn[j].lower_bound, xp[j].upper_bound, xn[j].upper_bound = 0, 0, inf, inf
			constraint[j]._constraint_sense = 'E'

		for i in metU:
			# S_ij x m_ie for each i
			m[i].add_metabolites({constraint[j]: j._metabolites[i] for j in list(i._reaction) if j in rxnK})
			m[i].upper_bound = inf

		#Add constraint to fix the total inconsistency for each element
		constraint_minIncon = Metabolite('minIncon')
		constraint_minIncon._bound = inf
		constraint_minIncon._constraint_sense = 'L'
		for j in rxnK:
			xp[j].add_metabolites({constraint_minIncon: 1})
			xn[j].add_metabolites({constraint_minIncon: 1})

		met_model.add_reactions([m[i] for i in metU])
		met_model.add_reactions([xp[j] for j in rxnK])
		met_model.add_reactions([xn[j] for j in rxnK])

		met_model.constraints[constraint_minIncon.id].lb = neg_inf #to avoid error

		objective_dict_minIncon = {xp[j]: 1 for j in rxnK}
		objective_dict_minIncon.update({xn[j]: 1 for j in rxnK})


		for e in ele:
			ct += 1
			print("Optimizing for %d / %d element ... %s" %(ct, len(ele), datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
			# metModelJ = Model('min/max ' + e)
			# metModelJ.solver = self.__solver
			infeasJ, boundJ, objJ, rhsJ = {}, {}, {}, {}
			# constraint = {j: Metabolite(j.id + ',' + e) for j in rxnK}
			# m = {i: Reaction('m_' + i.id + ',' + e) for i in metU}
			# xp = {j: Reaction('xp_' + j.id + ',' + e) for j in rxnK}
			# xn = {j: Reaction('xn_' + j.id + ',' + e) for j in rxnK}
			# for j in rxnK:
				#RHS for each constraint: -sum(S_ij * m^known_ie) (obsolete)
				# constraint[j]._bound = -sum([S_ij * metK[i].elements[e] for i, S_ij in j._metabolites.items() if i in metK and e in metK[i].elements])
				# constraint[j]._constraint_sense = 'E'
				#x_pos - x_neg for each constraint
				# xp[j].add_metabolites({constraint[j]: 1})
				# xn[j].add_metabolites({constraint[j]: -1})
				# xp[j].lower_bound, xn[j].lower_bound, xp[j].upper_bound, xn[j].upper_bound = 0, 0, inf, inf
			for i in metU:
				# # S_ij x m_ie for each i
				# m[i].add_metabolites({constraint[j]: j._metabolites[i] for j in list(i._reaction) if j in rxnK})
				# m[i].upper_bound = inf
				m[i].lower_bound = neg_inf if e == 'Charge' else 0
			#add reactions into the model
			# metModelJ.add_reactions([m[i] for i in metU])
			# metModelJ.add_reactions([xp[j] for j in rxnK])
			# metModelJ.add_reactions([xn[j] for j in rxnK])
			#set the objective function
			# objective_dict = {xp[j]: 1 for j in rxnK}
			# objective_dict.update({xn[j]: 1 for j in rxnK})
			met_model.constraints[constraint_minIncon.id].ub = inf
			set_objective(met_model, objective_dict_minIncon)
			#set RHS: -sum(S_ij * m^known_ie)
			for j in rxnK:
				met_model.constraints[constraint[j].id].ub = inf #ub must be set to be > lb to avoid error
				rhsJ[j] = -sum([S_ij * metK[i].elements[e] for i, S_ij in j._metabolites.items() if i in metK and e in metK[i].elements])
				met_model.constraints[constraint[j].id].lb, met_model.constraints[constraint[j].id].ub = rhsJ[j], rhsJ[j]
				constraint[j]._bound = rhsJ[j]

			#Solve for minimum inconsistency
			kwargs['objective_sense'] = 'minimize'
			sol = met_model.optimize(**kwargs)
			solStatJ = 'minIncon'
			infeasJ[solStatJ] = solution_infeasibility(met_model, sol)
			if sol.fluxes is None:
				boundJ[solStatJ] = float('nan')
			else:
				boundJ[solStatJ] = sum([sol.fluxes[j.id] for j in chain(xp.values(), xn.values())])
			if not infeasJ[solStatJ] <= feasTol:
				#infeasible (should not happen)
				infeasJ['minMw'], infeasJ['maxMw'] = inf, inf
				solStatJ = 'infeasible'
				objJ['minIncon'], objJ['minMw'], objJ['maxMw'] = (float('nan') for i in range(3))
				for s in ['minIncon','minMw','maxMw']:
					solution[s]['m'][e] = {i: float('nan') for i in metU}
					solution[s]['xp'][e] = {j: float('nan') for j in rxnK}
					solution[s]['xn'][e] = {j: float('nan') for j in rxnK}
			else:
				#Feasible. Store the solution
				solution[solStatJ]['m'][e] = {i: sol.fluxes[m[i].id] for i in metU}
				solution[solStatJ]['xp'][e] = {j: sol.fluxes[xp[j].id] for j in rxnK}
				solution[solStatJ]['xn'][e] = {j: sol.fluxes[xn[j].id] for j in rxnK}
				objJ[solStatJ] = sol.objective_value
				# #Add constraint to fix the total inconsistency for each element
				# constraint_minIncon = Metabolite('minIncon_'+e)
				# constraint_minIncon._bound = round(boundJ['minIncon'], digitRounded)
				# constraint_minIncon._constraint_sense = 'L'
				# for j in rxnK:
				# 	xp[j].add_metabolites({constraint_minIncon: 1})
				# 	xn[j].add_metabolites({constraint_minIncon: 1})

				# met_model.constraints[constraint_minIncon.id].lb = neg_inf #to avoid error
				met_model.constraints[constraint_minIncon.id].ub = round(boundJ['minIncon'], digitRounded)
				#reset the objective function to minimize molecular weight
				objective_dict = {m[metI]: Formula(formula=e).mw}
				set_objective(met_model, objective_dict)
				solStatJ = 'minMw'
				kwargs['objective_sense'] = 'minimize'
				eps0 = 1e-6
				while True:
					sol = met_model.optimize(**kwargs)
					infeasJ[solStatJ] = solution_infeasibility(met_model, sol)
					if infeasJ[solStatJ] <= feasTol or eps0 > 1e-4 + 1e-8:
						break
					eps0 *= 10
					#rounding to avoid infeasibility due to numerical issues
					met_model.constraints[constraint_minIncon.id].ub = round(boundJ['minIncon'] * (1 + eps0), digitRounded)

				boundJ[solStatJ] = eps0
				if infeasJ[solStatJ] <= feasTol:
					#Feasible. Store the solution
					solution[solStatJ]['m'].update({e: {i: sol.fluxes[m[i].id] for i in metU}})
					solution[solStatJ]['xp'].update({e: {j: sol.fluxes[xp[j].id] for j in rxnK}})
					solution[solStatJ]['xn'].update({e: {j: sol.fluxes[xn[j].id] for j in rxnK}})
					objJ[solStatJ] = sol.objective_value
				else:
					#infeasible, should not happen
					objJ[solStatJ] = float('nan')
					solution['minMw']['m'][e] = {i: float('nan') for i in metU}
					solution['minMw']['xp'][e] = {j: float('nan') for j in rxnK}
					solution['minMw']['xn'][e] = {j: float('nan') for j in rxnK}

				#maximize molecular weight
				solStatJ = 'maxMw'
				#reset the bound for total inconsistency
				met_model.constraints[constraint_minIncon.id].ub = round(boundJ['minIncon'], digitRounded)
				kwargs['objective_sense'] = 'maximize'
				eps0 = 1e-6
				while True:
					sol = met_model.optimize(**kwargs)
					infeasJ[solStatJ] = solution_infeasibility(met_model, sol)
					if infeasJ[solStatJ] <= feasTol or eps0 > 1e-4 + 1e-8:
						break
					eps0 *= 10
					#rounding to avoid infeasibility due to numerical issues
					met_model.constraints[constraint_minIncon.id].ub = round(boundJ['minIncon'] * (1 + eps0), digitRounded)

				boundJ[solStatJ] = eps0
				if infeasJ[solStatJ] <= feasTol:
					#Feasible. Store the solution
					solution[solStatJ]['m'].update({e: {i: sol.fluxes[m[i].id] for i in metU}})
					solution[solStatJ]['xp'].update({e: {j: sol.fluxes[xp[j].id] for j in rxnK}})
					solution[solStatJ]['xn'].update({e: {j: sol.fluxes[xn[j].id] for j in rxnK}})
					objJ[solStatJ] = sol.objective_value
				else:
					#infeasible, should not happen
					objJ[solStatJ] = float('nan')
					solution['maxMw']['m'][e] = {i: float('nan') for i in metU}
					solution['maxMw']['xp'][e] = {j: float('nan') for j in rxnK}
					solution['maxMw']['xn'][e] = {j: float('nan') for j in rxnK}

			#store data
			infeas[e] = infeasJ
			bound[e] = boundJ
			obj[e] = objJ
			# met_model[e] = met_model
			sol_stat[e] = solStatJ
			rhs[e] = rhsJ

		#summarize the final solution state
		if any([k == 'infeasible' for k in sol_stat.values()]):
			print('Failure: no feasible solution can be found.')
			solFinal = 'infeasible'
		else:
			if all([obj[e]['minMw'] == obj[e]['minMw'] for e in ele]):
				if all([obj[e]['maxMw'] == obj[e]['maxMw'] for e in ele]):
					solFinal = 'minMw + maxMw'
				else:
					solFinal = 'minMw only'
			else:
				if all([obj[e]['maxMw'] == obj[e]['maxMw'] for e in ele]):
					solFinal = 'maxMw only'
				else:
					solFinal = 'minIncon'

		mip_info = min_incon_parsi_info()
		mip_info.infeas, mip_info.bound, mip_info.obj, mip_info.solution, \
		mip_info.met_model, mip_info.final, mip_info.sol_stat \
		= infeas, bound, obj, solution, met_model, solFinal, sol_stat
		mip_info.rhs = rhs
		mip_info.mw_range = (sum([obj[e]['minMw'] for e in ele]), sum([obj[e]['maxMw'] for e in ele]))
		mip_info.formulae = tuple([formula_dict2str({e: solution[s]['m'][e][metI] for e in ele}) for s in ['minMw', 'maxMw']])
		return mip_info


	@property
	def solver(self):
		'''solver used. Default to be model.solver of the input cobra model.
		'''
		return self.__solver

	@solver.setter
	def solver(self, choice):
		'''change the optimization solver
		'''
		self.model.solver = choice
		self.__solver = self.model.solver
