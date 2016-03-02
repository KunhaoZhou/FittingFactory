import py_fit
import matplotlib.pyplot as plt
import numpy as np
from pylab import *

# register before use
py_fit.register_fit_factory()

# random x
x = 4.5 * np.random.randn(1000) 
x.sort()

# gaussian
yi = 10* np.exp(-0.5*((x-2)/5)**2)-6; 
[params, y_vals]= py_fit.fit(x.tolist(),yi.tolist(), "gaussian", "amoeba")
print("fitted gaussian model parameters: \n", params)
plt.plot(x,yi,"*",ms=10)
plt.plot(x,y_vals)
show()

# exponential
yi=2*np.exp(0.5*x)+10 
[params,y_vals]=py_fit.fit(x.tolist(),yi.tolist(), "exponential", "levenberg_marquardt");
print("fitted exponential model parameters: \n", params)
plt.plot(x,yi,"*",ms=10)
plt.plot(x,y_vals)
show()

# polynomial
yi=1*(x)+2*(x*x)+10*(x*x*x)+np.random.randn(1)
[params,y_vals]=py_fit.fit(x.tolist(),yi.tolist(), "polynomial", "svd");
print("fitted polynomial model parameters: \n", params)
plt.plot(x,yi,"*",ms=10)
plt.plot(x,y_vals)
show()

# radial basis
yi=1*np.exp(-(x -1)** 2 / 1)-2*np.exp(-(x -2)** 2 / 1)-3*np.exp(-(x -3)** 2 / 1)-3
[params,y_vals] = py_fit.fit(x.tolist(),yi.tolist(), "radial_basis", "levenberg_marquardt"); 
print("fitted radial_basis model parameters: \n", params)
plt.plot(x,yi,"*",ms = 10)
plt.plot(x,y_vals)
show()

# remove after use
py_fit.remove_fit_factory()
