import math as m
import numpy as np
import matplotlib.pyplot as pt

a = 3 #cm

da = 0.25 #cm

K = 1e-3

dt = 0.100

N = (int) (3*a/da) + 1

############### Initializing Variables ################

X = [[0 for i in range(0,N)]for n in range(0,N)]
Y = [[0 for i in range(0,N)]for n in range(0,N)]
T = [[0 for i in range(0,N)]for n in range(0,N)]
                                           
i = 0
j = 0
for i in range (0,N):
    for j in range (0,N):
        X[i][j] = da*j
        Y[i][j] = da*i

############ Creating Coefficient Matrix ####################

NC = N*N
R = K*dt/(2*da*da)

CM = [[0 for i in range(0,NC)]for n in range(0,NC)]
for j in range(0,N):
    for i in range(0,N):
        ii = j*N + i 
        if (X[i][j] > a) & (X[i][j] < 2*a) & (X[i][j] != 0) & (X[i][j] != (3*a)) & (Y[i][j] != 0) & (Y[i][j] != 3*a):
            CM[ii][ii] = 1 + 2*R + 2*R
            CM[ii][ii-1] = -R
            CM[ii][ii+1] = -R
            CM[ii][ii-N] = -R
            CM[ii][ii+N] = -R
        elif (Y[i][j] > a) & (Y[i][j] < 2*a) & (Y[i][j] != 0) & (Y[i][j] != 3*a)  & (X[i][j] != 0) & (X[i][j] != 3*a):
            CM[ii][ii] = 1 + 2*R + 2*R
            CM[ii][ii-1] = -R
            CM[ii][ii+1] = -R
            CM[ii][ii-N] = -R
            CM[ii][ii+N] = -R
            
        ################ Bottom
        
        elif (Y[i][j] > a) & (Y[i][j] < 2*a) & (X[i][j] == 0):
            CM[ii][ii] = 1 + 2*R + 2*R
            CM[ii][ii-1] = -R
            CM[ii][ii+1] = -R
            CM[ii][ii+N] = -2*R
        
        ################ Top
        
        elif (Y[i][j] > a) & (Y[i][j] < 2*a) & (X[i][j] == 3*a):
            CM[ii][ii] = 1 + 2*R + 2*R
            CM[ii][ii-1] = -R
            CM[ii][ii+1] = -R
            CM[ii][ii-N] = -2*R
        
        #################### Left
        
        elif (X[i][j] > a) & (X[i][j] < 2*a) & (Y[i][j] == 0):
            CM[ii][ii] = 1 + 2*R + 2*R
            CM[ii][ii+1] = -2*R
            CM[ii][ii-N] = -R
            CM[ii][ii+N] = -R
            
        #################### Right
        
        elif (X[i][j] > a) & (X[i][j] < 2*a) & (Y[i][j] == (3*a)):
            CM[ii][ii] = 1 + 2*R + 2*R
            CM[ii][ii-1] = -2*R
            CM[ii][ii-N] = -R
            CM[ii][ii+N] = -R
            
Time = 0

def ApplyBC(Time):
    for j in range(0,N):
        for i in range(0,N):
            ###### X direction
            if (X[i][j] == a) & (Y[i][j] <= a):
                T[i][j] = m.sin(Time)
            elif (X[i][j] == 2*a) & (Y[i][j] <= a):
                T[i][j] = m.cos(Time)
            elif (X[i][j] == 2*a) & (Y[i][j] >= 2*a):
                T[i][j] = m.sin(Time)
            elif (X[i][j] == a) & (Y[i][j] >= 2*a):
                T[i][j] = m.cos(Time)
                
            ###### Y Direction
            if (Y[i][j] == a) & (X[i][j] <= a):
                T[i][j] = m.sin(Time)
            elif (Y[i][j] == 2*a) & (X[i][j] <= a):
                T[i][j] = m.cos(Time)
            elif (Y[i][j] == 2*a) & (X[i][j] >= 2*a):
                T[i][j] = m.sin(Time)
            elif (Y[i][j] == a) & (X[i][j] >= 2*a):
                T[i][j] = m.cos(Time)
                
SM = [0 for i in range(0,NC)]

def ApplySolM():
    
    for j in range(0,N):
        for i in range(0,N):
            ii = j*N + i
            if (X[i][j] > a) & (X[i][j] < 2*a) & (X[i][j] != 0) & (X[i][j] != (3*a)) & (Y[i][j] != 0) & (Y[i][j] != 3*a):
                SM[ii] = T[i][j]*(1-4*R) +T[i-1][j]*R +T[i+1][j]*R+T[i][j-1]*R+T[i][j+1]*R
                    
            elif (Y[i][j] > a) & (Y[i][j] < 2*a) & (Y[i][j] != 0) & (Y[i][j] != 3*a)  & (X[i][j] != 0) & (X[i][j] != 3*a):
                SM[ii] = T[i][j]*(1-4*R) +T[i-1][j]*R +T[i+1][j]*R+T[i][j-1]*R+T[i][j+1]*R
                    
ApplyBC(Time)    
colorinterpolation = 50
#colourMap = pt.cm.jet # try: pt.cm.coolwarm
colourMap = pt.cm.coolwarm
pic = 0
pt.figure(1)
pt.title("Time = 0")
pt.contourf(X, Y, T,colorinterpolation)
pt.colorbar() # set colorbar
pt.show()
while (Time <= 4):
    Time += dt
    Temp = []
    ApplyBC(Time)
    for i in range(0,N):
        Temp = np.append(Temp,T[i][:])
    SM = Temp
    ApplySolM()
    maxres = 2e-3
    while (maxres < 1e-3):
        
        for i in range(0,NC):
            Temperature = SM[i]
            TT = Temp[i]
            for j in range(0,NC):
                if(i != j):
                    Temperature = Temperature + CM[i][j]*Temp[j]
            Temp[i] = Temperature/CM[i][i]
            Temp[i] = TT + 1.95*(Temp[i] - TT)
            res = abs(Temp[i] - TT)
            if (res>maxres):
                maxres = res
    p = 0
    for i in range(0,NC,N):
        T[p][:] = Temp[i:(i+N)]
        p += 1
    if (Time == 0.6):
        pt.figure(2)
        pt.title("Time = 0.6")
        pt.contourf(X, Y, T,colorinterpolation)
        pt.colorbar() # set colorbar
        pt.show()
    if (Time >= 3.1) & (pic == 0):
        pt.figure(3)
        pt.title("Time = 3.14159")
        pt.contourf(X, Y, T,colorinterpolation)
        pt.colorbar() # set colorbar
        pt.show()
        pic = 1

        
            
