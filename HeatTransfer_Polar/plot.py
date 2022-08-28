##############
# Author: Neelay Doshi
##############

##############
# Clearing variables 
#from IPython import get_ipython
#get_ipython().magic('reset -sf')

import numpy as np
import matplotlib.pyplot as pt

##-- Generate Data -----------------------------------------
## Using linspace so that the endpoint of 360 is included...
#azimuths = np.radians(np.linspace(0, 360, 20))
#zeniths = np.arange(0, 70, 10)
#
#r, theta = np.meshgrid(zeniths, azimuths)
#values = np.random.random((azimuths.size, zeniths.size))
#
##-- Plot... ------------------------------------------------
#fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
#ax.contourf(theta, r, values)
#
#plt.show()



out     = np.loadtxt('Dets.txt', dtype=int)
T_new   = np.loadtxt('Temp.txt', dtype=float)

#out     = np.loadtxt('Results/out_a_10_10.txt', dtype=int)
#T_new   = np.loadtxt('Results/T_new_a_10_10.txt', dtype=float)

#out     = np.loadtxt('Results/out_a_5_5.txt', dtype=int)
#T_new   = np.loadtxt('Results/T_new_a_5_5.txt', dtype=float)

#out     = np.loadtxt('Results/out_b_4_4.txt', dtype=int)
#T_new   = np.loadtxt('Results/T_new_b_4_4.txt', dtype=float)

problem = out[0]
n_phi   = out[1]
n_r     = out[2]
T_origin= T_new[0] #out[2]
Time    = 100


#-- Generate Data -----------------------------------------
# Using linspace so that the endpoint of 360 is included...
azimuths = np.radians(np.linspace(0, 360, n_phi))
zeniths = np.linspace(0, 5e-2, n_r)

r, theta    = np.meshgrid(zeniths, azimuths)
values      = np.zeros( ( len(azimuths), len(zeniths) ) )

count = 1
for i in range(len(azimuths)-1):
    values[i][0] = T_origin
    for j in range(len(zeniths)):
        
        if (j != 0):
            values[i][j] = T_new[count]
            count += 1

values[-1][:] = values[0][:]
#print(values)

#-- Plot... ------------------------------------------------
colourMap = pt.cm.seismic
fig, ax = pt.subplots(subplot_kw=dict(projection='polar'))

strg = 'A';
if problem == 2:
    strg = 'B'

pt.contourf(theta, r, values, 50, cmap= 'coolwarm')
pt.colorbar() # set colorbar
pt.title("Contour of Temperature")
#pt.show()
string = 'Plots/Countour_' + strg + str(n_phi) + 'x' + str(n_r) + '.png'
pt.savefig(string)




