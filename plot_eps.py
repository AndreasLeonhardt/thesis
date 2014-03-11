import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import array
import os;

filename = sys.argv[1];
A=loadtxt(filename, skiprows=0);

[nx,ny]=A.shape

X=[np.pi*(-1.+2.*n/nx) for n in range(nx)]
Y=[np.pi*(-1.+2.*n/ny) for n in range(ny)]


[X,Y]=np.meshgrid(X,Y)

f = figure();
ad = f.add_subplot(111,projection='3d');

ad.plot_wireframe(X,Y,A)
#ad.plot_wireframe(X,Y,Ep,color='red')

#xlabel('j');
#ylabel('E');
#title('blocking results');




show(block=True)

quit()
