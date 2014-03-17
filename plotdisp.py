from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import array
import os;

filename = sys.argv[1];
A=loadtxt(filename, skiprows=0);

nx=int(A[-1,0])
ny=int(A[-1,0])

X=zeros((nx,ny));
Y=zeros((nx,ny));
Em=zeros((nx,ny));
Ep=zeros((nx,ny))
eps=zeros((nx,ny))

for i in range(nx):
	for j in range(ny):
		X[i,j]=(A[i*(ny+1),2]);
		Y[i,j]=(A[j,3]);
		Em[i,j]=(A[i+j*(nx+1),5]);
		Ep[i,j]=(A[i+j*(nx+1),6]);
		eps[i,j]=A[i+j*(nx+1),4]





f = figure();

ad = f.add_subplot(111,projection='3d');
ad.plot_wireframe(X,Y,eps)
#ad.plot_wireframe(X,Y,Ep,color='red')

#ad = f.add_subplot(111)
#p=ad.contour(X,Y,eps)
#cb=f.colorbar(p,ax=ad)

#xlabel('j');
#ylabel('E');
#title('blocking results');




show(block=True)

quit()
