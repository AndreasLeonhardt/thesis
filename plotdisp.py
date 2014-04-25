import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import array
import os

filename = sys.argv[1]
A=loadtxt(filename, skiprows=0)

nx=int(A[-1,0])
ny=int(A[-1,0])

X=np.zeros((nx+1,ny+1))
Y=np.zeros((nx+1,ny+1))
Em=np.zeros((nx+1,ny+1))
Ep=np.zeros((nx+1,ny+1))
eps=np.zeros((nx+1,ny+1))

for i in range(nx):
	for j in range(ny):
		X[i,j]=(A[i*(ny+1),2])/(2*np.pi)
		Y[i,j]=(A[j,3])/(2*np.pi)
		Em[i,j]=(A[i+j*(nx+1),5])
		Ep[i,j]=(A[i+j*(nx+1),6])
		eps[i,j]=A[i+j*(nx+1),4]#*.259#eV

for j in range(ny):
	X[nx,j]=-X[0,j]
	eps[nx,j]=eps[0,j]
	Y[nx,j]=Y[0,j]

for i in range(nx+1):
	X[i,ny]=X[i,0];
	Y[i,ny]=-Y[i,0];
	eps[i,ny]=eps[i,0];





t=np.arange(np.floor(eps.min()),np.ceil(eps.max()),1.)

f = figure()

#ad = f.add_subplot(111,projection='3d')
#ad.plot_wireframe(X,Y,eps)
#ad.plot_wireframe(X,Y,Ep,color='red')

ad = f.add_subplot(111)
ad.arrow(-.405,-.095,.5-0.04/np.sqrt(2),.5-0.04/np.sqrt(2),head_width=.02,head_length=.04,fc="k", ec="k" )
p=ad.contour(X,Y,eps,t)
cb=f.colorbar(p,ax=ad)
plt.xticks(np.arange(-.5,.6,.25))
plt.yticks(np.arange(-.5,.6,.25))
#xlabel('j')
#ylabel('E')
#title('blocking results')




show(block=True)

quit()
