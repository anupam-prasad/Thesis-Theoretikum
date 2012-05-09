import pickle as pickle
import numpy as np
import matplotlib.pyplot as plt

f1=open('test5results3','r')
f2=open('test5results3_unmodified','r')

y1=pickle.load(f1)
y2=pickle.load(f2)

plotval1=abs(y1[3:,:])
plotval2=abs(y2[3:,:])

n=np.linspace(0,14,15);

for k in range(4):
	mom_index=k+4
	axlim1=min(plotval1[:,mom_index])
	axlim2=1.01*max(plotval1[:, mom_index])
	yexact=abs(y1[0,mom_index]) * np.ones(15)
	plt.subplot(2,2,k+1)
	plt.plot(n,plotval1[:,mom_index],'ro',n,yexact,'b--')
	plt.axis([-1, 15, axlim1, axlim2])

plt.show()


#print plotval2[:,0]
#exit('here')
for k in range(4):
	mom_index=k+4
	axlim1=min(plotval2[:,mom_index])
	axlim2=max(plotval2[:,mom_index])
	yexact=abs(y2[0,mom_index]) * np.ones(15)
	plt.subplot(2,2,k+1)
#	plt.plot(n[:4],np.log10(plotval2[:4,mom_index]),'ro',n[:4],np.log10(yexact),'b--')
#	plt.axis([-1, 15, axlim1, axlim2])
	plt.plot(n,plotval2[:,mom_index],'ro',n,yexact,'b--')
	plt.axis([-1, 15, axlim1, axlim2])

plt.show()

	




