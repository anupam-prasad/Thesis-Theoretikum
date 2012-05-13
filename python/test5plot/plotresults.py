import pickle as pickle
import numpy as np
import matplotlib.pyplot as plt

f1=open('test5results_new','r')
f2=open('test5results_new_unmod','r')

y1=pickle.load(f1)
y2=pickle.load(f2)


plotval1=abs(y1[3:,:])
plotval2=abs(y2[3:,:])

n=np.linspace(0,19,20);

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

for k in range(3):
	mom_index=k

	axlim1=.9*min(plotval1[:,mom_index])
	axlim2=1.01*max(plotval1[:, mom_index])
	yexact=abs(y1[0,mom_index]) * np.ones(20)
	plt.subplot(1,2,1)
	plt.plot(n,plotval1[:,mom_index],'ro',n,yexact,'b--')
	plt.plot(n,plotval1[:,mom_index],'r--')
#	plt.axis([-1, 15, axlim1, axlim2])

	axlim1=min(plotval2[:,mom_index])
	axlim2=max(plotval2[:,mom_index])
	yexact=abs(y2[0,mom_index]) * np.ones(20)
	plt.subplot(1,2,2)
	plt.plot(n,plotval2[:,mom_index],'ro',n,yexact,'b--')
	plt.plot(n,plotval2[:,mom_index],'r--')
#	plt.axis([-1, 15, axlim1, axlim2])
	plt.suptitle('E_n=%s*pi^2/200'%(mom_index*mom_index))
	plt.figtext(.45,.02,'x-axis = no. of iterations')
	plt.show()

	




