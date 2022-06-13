# %%
import numpy as np
import matplotlib.pyplot as plt

# b=loadtxt("evs.dat")
# c=loadtxt("evcs.dat")

i = 2

if i == 0:
    a = np.loadtxt("../build/dataBaseFlow.txt")
    plt.plot(a[:, 0], "-x")
    plt.plot(a[:, 1], "-x")
    plt.plot(a[:, 2], "-x")
    plt.show()

if i == 1:
    a = np.loadtxt("../build/dataBaseFlowPress.txt")
    plt.plot(a[:, 0], "-x")
    plt.plot(a[:, 1], "-x")
    plt.plot(a[:, 2], "-x")
    plt.plot(a[:, 3], "-x")
    plt.plot(a[:, 4], "-x")
    plt.show()

if i == 2:
    a = np.loadtxt("../build/evs.dat")
    plt.plot(a[:,0],a[:,1], "x")
    plt.grid()
    plt.show()
    print(np.max(a[:,1]))
    print(a[(np.argmax(a[:,1])),0])

if i == 3:
    a = np.loadtxt("../build/evcs.dat")
    plt.plot(a[0,4:-1:4], "-x")
    plt.plot(a[1,4:-1:4], "-x")
    plt.show()
# %%
