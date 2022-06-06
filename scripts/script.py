#%%
import numpy as np
import matplotlib.pyplot as plt

#b=loadtxt("evs.dat")
#c=loadtxt("evcs.dat")

i=1

if i==0:
	a=np.loadtxt("../build/dataBaseFlow.txt")
	plt.plot(a[:,0],"-x")
	plt.plot(a[:,1],"-x")
	plt.plot(a[:,2],"-x")
	plt.show()

if i==1:
	a=np.loadtxt("../build/dataBaseFlowPress.txt")
	plt.plot(a[:,0],"-x")
	plt.plot(a[:,1],"-x")
	plt.plot(a[:,2],"-x")
	plt.plot(a[:,3],"-x")
	plt.plot(a[:,4],"-x")
	plt.show()
# %%
