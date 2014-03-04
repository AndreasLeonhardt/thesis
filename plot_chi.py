# python script for plotting results from Susceptibility
# call with data file as argument
import  numpy
import pylab
from subprocess import call
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import array
import os

# load data from file, using the first argument as filename.
filename = pylab.sys.argv[1]
A=pylab.loadtxt(filename, skiprows=2)

# get max coordinates (q,w)
# (might need modification for offset, e.g. not starting with w=0.)
tmax=int(A[-1,0])+1
n_w =int(len(A[:,0])/(tmax))

# specifier for type of plot, as written from Susceptibility.
x=4
xb=5
y=6
yb=7
z1=8
z1b=9
z2=10
z2b=11


# values for old files, 
# e.g. where yb=y, z1b=z1, z2b=z2
#uncommment the following line
#yb=6; z1=7; z1b=7; z2=8; z2b=8



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
T=numpy.array(range(tmax)).astype(float)
W=numpy.array(range(n_w )).astype(float)
Qx=[]
Qy=[]

# fill with corresponding coordinates/values
for t in range (0,tmax):
	T[t]=A[t*n_w,0]
# plotted values are selected (x,xb,y,z1,z2) /calculated (chi, denominator of chi) here.
	for w in range(0,n_w):
		tw=t*n_w+w
		
		# choose x,xb,y,z1 or z2
		#Val[t,w]=A[tw,z2]
		
		# calculate denominator of chi
		#Val[t,w]=(1.0+lmb*(A[tw, x]+A[tw, y]))*(1.0+lmb*(A[tw,xb]+A[tw,yb]))-lmb**2*(A[tw,z1]+A[tw,z2])*(A[tw,z1b]+A[tw,z2b])
		
		# calculate chi
		
		Val[t,w]=((-(A[tw, x]+A[tw, y])*(1.0+lmb*(A[tw,xb]+A[tw,yb]))+lmb*(A[tw,z1]+A[tw,z2])*(A[tw,z1b]+A[tw,z2b]))
		/((1.0+lmb*(A[tw, x]+A[tw, y]))*(1.0+lmb*(A[tw,xb]+A[tw,yb]))-lmb**2*(A[tw,z1]+A[tw,z2])*(A[tw,z1b]+A[tw,z2b])))
		
		# define X,Xb,Y,Z1u,Z2u,Z1d,Z2d for usage below, spin is differentiated by u and d only when it matters.
		K = +A[tw, x] -A[tw, y]
		Kb= +A[tw,xb] -A[tw,yb]
		N = +A[tw,z1] -A[tw,z2] 
		Nb= +A[tw,z1b] -A[tw,z2b]

		# chi^zz (denominator fist, then the whole expression, using the denominator (comment carefully))
		#Val[t,w]= ((1.0+lmb**2*(-K**2+N**2))*(1.0+lmb**2*(-Kb**2+N**2))+lmb**4*(K*N-N*Kb))
		# calculate enumerator, dividing by the previous expression. 
		#Val[t,w]=2*((1+lmb**2*(N**2-Kb**2))*(K-lmb*K**2)-lmb**2*(K-Kb)*(N**2-lmb*K*N**2)-lmb**3*K*(K-Kb)*N**2-lmb*N**2*(1+lmb**2*(N**2-K**2)))/Val[t,w]
	
	Qx.append(A[t*n_w,1])
	Qy.append(A[t*n_w,2])


for w in range(n_w):
	W[w]=float(A[w,3])
    
# create 2D arrays for plot function using meshgrid (see doc)
W,T=numpy.meshgrid(W,T)


# create figure, subfigure and plot finally
fig = plt.figure()
#ad = fig.add_subplot(1,1,1,projection='3d')
ad = fig.add_subplot(1,1,1)

#p= ad.plot_wireframe(T,W,Val)
p=ad.imshow(Val.transpose(),extent = [T.min(),T.max(),W.min(),W.max()], aspect='auto',origin = 'lower',cmap = 'RdBu')
p.set_clim(-100000,100000)
#p= ad.contour(T,W,Val)
cb = fig.colorbar(p,ax=ad)

# limit z axis, since results may diverge
#ad.set_zlim3d(-5000,5000)

# plot path through Brillouin zone
#ad= fig.add_subplot(1,1,1)
#p=ad.plot(Qx,Qy)
#ad.set_xlim(-numpy.pi,numpy.pi)
#ad.set_ylim(-numpy.pi,numpy.pi)

pylab.show(block=True)

quit()
