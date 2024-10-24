import numpy as np
from sympy import Symbol
from sympy import lambdify
from sympy import *
import matplotlib.pyplot as plt
#Nodes that want to be interpolated are given
x=np.array([0,0.6,0.9])
y=np.array([1,0.83,0.62])
lag_coef=[]
#Symbolic polinom is defined
def syml(t):
	lag=1.
	n=Symbol("n")
	for i in range(len(x)):
		if i!=t:
			lag*=(n-i)/(x[t]-x[i])
		else:
			continue
	return lag
print(lag_coef)
def symp():
	pol=1.
	for i in range(len(x)):
		pol+=syml(i)*y[i]
	return pol
print(factor(symp()))
#t=np.linspace(x[0],x[len(x)-1],100)
#plt.plot(t,pol(t),color="black",label="Interpolation")
#plt.scatter(x,y,color="red",label="Nodes")
#plt.ylim(0,2)
#plt.title("Newtons interpolation method")
#plt.legend()
#plt.grid()
#plt.show()
