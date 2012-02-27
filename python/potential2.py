import numpy as np

import scipy

def potential(x,potential_type,lb=0.,ub=10.):
	V=np.zeros(len(x))

	if potential_type=='qho':
		for k in range(len(x)):
			V[k]=x[k]*x[k]/2

	elif potential_type=='fwell':
		for k in range(len(x)):
			if x[k] < 6. and x[k] > 4.:	V[k]=-400.
			else:	V[k]=0

	elif potential_type=='infwell':
		for k in range(len(x)):
			if x[k] < ub and x[k] > lb:	V[k]=0.
		
	
	elif potential_type=='gaussian':
		for k in range(len(x)):
			V[k]=-10*(np.exp(-((x[k]+1)*(x[k]+1))/.05))
	
	return V

def potential2(x):
	V=np.zeros(len(x))
	for k in range(len(x)):
#		V[k]=x[k]*x[k]/2		
		if x[k] < 3. and x[k] > 0.:
#				V[k]=-10
#                               V[k]=-15*np.exp(-x[k]*x[k]/.05)
                                V[k]=-1/(abs(x[k])+.1)
#                               V[k]=-(np.exp(-abs(x[k]))) - 1/(abs(x[k])+.1)
#  				V[k]=0
		else:
			V[k]=0

	return V


def potential3(x):
	V=np.zeros(len(x))
	for k in range(len(x)):
		if x[k] < -.5 and x[k] > -1.5:
				V[k]=-10
		elif x[k] < 1.5 and x[k] > .5:
				V[k]=-5
		else:	V[k]=0

	return V
