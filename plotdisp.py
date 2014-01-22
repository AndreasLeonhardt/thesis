from numpy import *
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import array
import os;

filename = sys.argv[1];
A=loadtxt(filename, skiprows=0);

begin=0;
n=50;
treshold=1*n*n

X=zeros((n,n));
Y=zeros((n,n));
Em=zeros((n,n));
Ep=zeros((n,n));


for i in range(0,n):
	for j in range (0,n):
		X[i,j]=(A[i+j*n+treshold,1]);
		Y[i,j]=(A[i+j*n+treshold,2]);
		Em[i,j]=(A[i+j*n+treshold,3]);
		Ep[i,j]=(A[i+j*n+treshold,4]);






f = figure();
ad = f.add_subplot(111,projection='3d');

ad.plot_wireframe(X,Y,Em)
ad.plot_wireframe(X,Y,Ep,color='red')

#xlabel('j');
#ylabel('E');
#title('blocking results');




show(block=True)

quit()
