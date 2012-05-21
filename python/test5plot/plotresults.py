import pickle as pickle
import numpy as np
import matplotlib.pyplot as plt

f1=open('test5results_fwell_mod','r')
f2=open('test5results_fwell_unmod','r')

y1=pickle.load(f1)
y2=pickle.load(f2)


plotval1=abs(y1[3:,:])
plotval2=abs(y2[3:,:])

N=len(plotval1[:,1])
N2=len(plotval2[:,1])

n=np.linspace(0,N-1,N);
n2=np.linspace(0,N2-1,N2);

#for k in range(4):
#	mom_index=k+4
#	axlim1=min(plotval1[:,mom_index])
#	axlim2=1.01*max(plotval1[:, mom_index])
#	yexact=abs(y1[0,mom_index]) * np.ones(15)
#	plt.subplot(2,2,k+1)
#	plt.plot(n,plotval1[:,mom_index],'ro',n,yexact,'b--')
#	plt.axis([-1, 15, axlim1, axlim2])
#
#plt.show()
#
#
##print plotval2[:,0]
##exit('here')
#for k in range(4):
#	mom_index=k+4
#	axlim1=min(plotval2[:,mom_index])
#	axlim2=max(plotval2[:,mom_index])
#	yexact=abs(y2[0,mom_index]) * np.ones(15)
#	plt.subplot(2,2,k+1)
#	plt.plot(n,plotval2[:,mom_index],'ro',n,yexact,'b--')
#	plt.axis([-1, 15, axlim1, axlim2])
#
#plt.show()

for k in range(len(plotval1[0,:])):
	mom_index=k

	axlim1=.9*min(plotval1[:,mom_index])
	axlim2=1.01*max(plotval1[:, mom_index])
	yexact=abs(y1[0,mom_index]) * np.ones(N)
	plt.subplot(1,2,1)
	plt.plot(n,plotval1[:,mom_index],'ro',n,yexact,'b--')
	plt.plot(n,plotval1[:,mom_index],'r--')
	plt.ylabel('G - matrix element')
#	plt.axis([-1, 15, axlim1, axlim2])

	axlim1=min(plotval2[:,mom_index])
	axlim2=max(plotval2[:,mom_index])
	yexact=abs(y2[1,mom_index]) * np.ones(N2)
	plt.subplot(1,2,2)
	plt.plot(n2,plotval2[:,mom_index],'ro',n2,yexact,'b--')
	plt.plot(n2,plotval2[:,mom_index],'r--')
#	plt.axis([-1, 15, axlim1, axlim2])
	plt.suptitle('E_n=%s*pi^2/200'%(mom_index*mom_index))
	plt.figtext(.45,.02,'x-axis = no. of iterations')
	plt.show()

	




