from __future__ import print_function

class data_object(object):
	'''data_object has a representation summarizing all public attributes
	'''
	def __repr__(self):
		print('Attribute\tClass\tValue/Length')
		for k in dir(self):
			if k[0] != '_':
				val = eval("self." + k)
				s = ''
				if isinstance(val, (int, long, float, complex)):
					s = str(val)
				elif isinstance(val, (list, dict, tuple, set)):
					s = str(len(val))
				elif isinstance(val, str):
					s = val
				print('%s\t%s\t%s\t at 0x%x' %(k, type(val), s, id(val)))
		print('')
		return "<%s at 0x%x>" % (self.__class__.__name__, id(self))

class preprocess_data(data_object):
	'''Preprocess data:
	met_known: {i: GenericFormula(i)} for all metabolite i in the input model
	met_unknown: list of all unknown metabolite objects in the input model
	met_fill: {'X1Y2Z3': GenericFormula('X1Y2Z3')} for all filling metabolite (e.g. with formula 'X1Y2Z3')
	rxn_known: list of reaction objects in the input model used for balancing
	ele: list of elements existing in the chemical formulae of met_known
	ele_connect: list of connected elements in tuples
	met_fill_connect: {ele_connent: met_fill}, dictionary associating each filling metabolite met_fill to the connected sets of elements ele_connect.
	'''
	def __init__(self):
		self.met_known, self.met_unknown, self.met_fill, self.rxn_known, \
		self.ele, self.ele_connect, self.met_fill_connect = (None for i in range(7))

class min_incon_parsi_info(data_object):
	'''min_incon_parsi_info object summarizing the results of solving MIP.
		formulae: result formulae
		mw_range: (min MW, max MW) (computeMetRange only)
		rhs: {e: {i: RHS[i] for all metabolite i}} the RHS value in the MIP problem for each element solved. (computeMetRange only)
		infeas: infeasibility of each solve
		bound: bound used for total inconsistency or the relaxation value eps0 for each solve
		obj: objective function value for each solve
		solution: solution values for each type of variables (m, xp, xn, Ap, An)
		met_model: the MIP problems solved for each element/connected set of elements as individual cobra models
		final: the final status of the solution
		sol_stat: solution status for each element/connected set of elements
		sol_constrain: 'minIncon' or 'minFill' indicating the solution used to constrain the minimal formulae problem. (computeMetForm only)
	'''
	def __init__(self):
		self.formulae, self.mw_range, self.rhs, self.infeas, self.bound, self.obj, self.solution, self.met_model, \
		self.final, self.sol_stat, self.sol_constrain = (None for i in range(11))

class conserved_moiety_info(data_object):
	'''conserved_moiety_info object containing the conserved moiety information:
		cm: list of conserved moieties, each being a metabolite-coefficient dictionary
		cm_generic: list of entries in cm that are generic (no known metabolites involved)
		cm_generic_dict: {fomrula: cm_generic} dictionary. formula is the GenericFormula object with the defaulted or inputted name as the formula for the conserved moiety.
	'''
	def __init__(self):
		self.cm, self.cm_generic, self.cm_generic_dict = (None for i in range(3))

class result_structure(data_object):
	'''ResultStructure object with the following attributes:
		formulae: the final formulae for unknown metabolites
		S_fill: a dictionary of stoichiometric coefficients for filling metabolites
		rxn_balance: {rxn:{e: balance}} dictionary for the elemental balance for each element e of each reaction rxn
		model: a copy of the input cobra model with updated chemical formulae
		mip_info: min_incon_parsi_info object summarizing the results of solving MIP.
		cm_info: conserved_moiety_info object containing the conserved moiety information.
	'''
	def __init__(self):
		self.formulae, self.S_fill, self.rxn_balance, self.model, self.mip_info, self.cm_info = (None for i in range(6))
