import numpy as np

import scipy

def potential(x,potential_type,params=[]):
	V=np.zeros(len(x))

	if np.size(params) > 0:	
		lb=params[0]
		ub=params[1]
		pot_strength=params[2]

	if potential_type=='qho':
		for k in range(len(x)):
			V[k]=x[k]*x[k]/2

	elif potential_type=='fwell':
		for k in range(len(x)):
			if x[k] < ub and x[k] > lb:	V[k]=-pot_strength
			else:	V[k]=0

	elif potential_type=='infwell':
		for k in range(len(x)):
			V[k]=0.
		
	
	elif potential_type=='gaussian':
		center=(lb+ub)/2.
		sigma=(ub-lb)/2.
		for k in range(len(x)):
			V[k]=-pot_strength*(np.exp(-((x[k]-center)*(x[k]-center))/(2.*sigma*sigma)))
	
	elif potential_type=='gaussiancutoff':
		center=(lb+ub)/2.
		sigma=(ub-lb)/2.
		for k in range(len(x)):
			if x[k] < ub and x[k] > lb:	V[k]=-pot_strength*(np.exp(-((x[k]-center)*(x[k]-center))/(2.*sigma*sigma)))
			else:	V[k]=0
	
	return V

