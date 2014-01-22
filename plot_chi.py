# python script for plotting results from Susceptibility
# call with data file as argument
import  numpy
from pylab import *
from subprocess import call
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import array
import os

# load data from file, using the first argument as filename.
filename = sys.argv[1]
A=loadtxt(filename, skiprows=2)

# get max coordinates (q,w)
# (might need modification for offset, e.g. not starting with w=0.)
tmax=int(A[-1,0])+1
n_w =int(len(A[:,0])/(tmax))

# specifier for type of plot, as written from Susceptibility.
x=4
xb=5
y=6
z1=7
z2=8

# get external parameters
f = open(filename, 'r')
paramstext = f.readline()
f.close()
params=paramstext.split()
# assign them to variables
U=float(params[1])
N=float(params[9])
lmb=U/N


# create numpy arrays for coordinates and values
Val=numpy.ndarray((tmax,n_w))          
T=numpy.array(range(tmax))
W=numpy.array(range(n_w ))
Qx=[]
Qy=[]

# fill with corresponding coordinates/values
for t in range (0,tmax):
    T[t]=A[t*n_w,0]
# plotted values are selected (x,xb,y,z1,z2) /calculated (chi, denominator of chi) here.
    for w in range(0,n_w):
        # choose x,xb,y,z1 or z2
	#Val[t,w]=A[t*n_w+w,z1]
	# calculate denominator of chi
	#Val[t,w]=(1+lmb*(A[t*n_w+w, x]+A[t*n_w+w, y]))*(1+lmb*(A[t*n_w+w,xb]*A[t*n_w+w, y]))-lmb**2*(A[t*n_w+w,z1]+A[t*n_w+w,z1])**2
	# calculate chi
	Val[t,w]=log((-(A[t*n_w+w, x]+A[t*n_w+w, y])*(1.0+lmb*(A[t*n_w+w,xb]+A[t*n_w+w, y]))+lmb*(A[t*n_w+w,z1]+A[t*n_w+w,z2])**2)
		/((1.0+lmb*(A[t*n_w+w, x]+A[t*n_w+w, y]))*(1.0+lmb*(A[t*n_w+w,xb]*A[t*n_w+w, y]))-lmb**2*(A[t*n_w+w,z1]+A[t*n_w+w,z1])**2))

    Qx.append(A[t*n_w,1])
    Qy.append(A[t*n_w,2])


for w in range(n_w):
    W[w]=A[w,3]
# create 2D arrays for plot function using meshgrid (see doc)
T,W=numpy.meshgrid(W,T)


# create figure, subfigure and plot finally
fig = plt.figure()
ad = fig.add_subplot(1,1,1,projection='3d')
p= ad.plot_surface(T,W,Val)
# limit z axis, since results may diverge
#ad.set_zlim3d(-500,500)

# plot path through Brillouin zone
ad= fig.add_subplot(3,3,3)
p=ad.plot(Qx,Qy)
ad.set_xlim(-numpy.pi,numpy.pi)
ad.set_ylim(-numpy.pi,numpy.pi)

show(block=True)

quit()
