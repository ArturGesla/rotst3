import numpy as np
import matplotlib.pyplot as plt

a=loadtxt("dataBaseFlow.txt")
b=loadtxt("evs.dat")
c=loadtxt("evcs.dat")

i=0

if i==0:
	plt.plot(a[:,0],"-x")
	plt.plot(a[:,1],"-x")
	plt.plot(a[:,2],"-x")
	plt.show()