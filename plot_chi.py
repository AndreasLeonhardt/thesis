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
# where b indicates a bar (q+Q) and r and i the real and imaginary parts respectiviely
n_xr=4
n_xi=5
n_xbr=6
n_xbi=7
n_yr=8
n_yi=9
n_ybr=10
n_ybi=11
n_z1r=12
n_z1i=13
n_z1br=14
n_z1bi=15
n_z2r=16
n_z2i=17
n_z2br=18
n_z2bi=19


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
		x=complex(A[tw,n_xr],A[tw,n_xi])
		xb=complex(A[tw,n_xbr],A[tw,n_xbi])
		y=complex(A[tw,n_yr],A[tw,n_yi])
		yb=complex(A[tw,n_ybr],A[tw,n_ybi])
		z1=complex(A[tw,n_z1r],A[tw,n_z1i])
		z1b=complex(A[tw,n_z1br],A[tw,n_z1bi])
		z2=complex(A[tw,n_z2r],A[tw,n_z2i])
		z2b=complex(A[tw,n_z2br],A[tw,n_z2bi])

		# calculate denominator of chi
		#Val[t,w]=(1.0+lmb*(x+y))*(1.0+lmb*(xb+yb))-lmb**2*(z1+z2)*(z1b+z2b)
		
		# calculate chi (imaginary part of the susceptibility
		
		# transversal (+-)
		Val[t,w]=((-(x+y)*(1.0+lmb*(xb+yb))+lmb*(z1+z2)*(z1b+z2b))/((1.0+lmb*(x+y))*(1.0+lmb*(xb+yb))-lmb**2*(z1+z2)*(z1b+z2b))).imag
		
		# longotudinal (zz)
		# define X,Xb,Y,Z1u,Z2u,Z1d,Z2d for usage below, spin is differentiated by u and d only when it matters.
		K = +x-y
		Kb= +xb-yb
		N = -z1+z2
		Nb= -z1b+z2b

		# calculate imaginary part of longitudinal retarted susceptibility. Factor of 2 because of spin (change up and down)
		#Val[t,w]=(-2*((1.0-lmb**2*(Kb*Kb-Nb*N))*(K+lmb*K*K)+lmb**3*(N*Nb*K*Kb-K*K*N*Nb)+lmb**2*N*Nb*(Kb-K)*(1.0+lmb*K)-lmb*(1.0+lmb**2*(N*Nb-K*K))*N*Nb)/((1.0+lmb**2*(N*Nb-Kb*Kb))*(1.0+lmb**2*(-K*K+N*Nb))-lmb**4*(-K*K*Nb*N+2*K*Kb*N*Nb-Kb*Kb*Nb*N))).imag
	
	Qx.append(A[t*n_w,1])
	Qy.append(A[t*n_w,2])


for w in range(n_w):
	# scale frequency to Heisenberg model predictions, that is with 4t^2/U 
	#(since w is in units of t or t=1, that would be divided by U/t, which is U here)
	W[w]=float(A[w,3])/4*U
    
# create 2D arrays for plot function using meshgrid (see doc)
W,T=numpy.meshgrid(W,T)


# create figure, subfigure and plot finally
fig = plt.figure()
#ad = fig.add_subplot(1,1,1,projection='3d')
ad = fig.add_subplot(1,1,1)

#p= ad.plot_wireframe(T,W,Val)
#p=ad.imshow(Val.transpose(),extent = [T.min(),T.max(),W.min(),W.max()], aspect='auto',origin = 'lower',cmap = 'RdBu') # RdBu, Blue
#p.set_clim(-5000,5000)
#p= ad.contour(T,W,Val)
#cb = fig.colorbar(p,ax=ad)

p=ad.plot(W[0,:],Val[0,:],W[96,:],Val[96,:])

# limit z axis, since results may diverge
#ad.set_zlim3d(-5000,5000)

# plot path through Brillouin zone
#ad= fig.add_subplot(1,1,1)
#p=ad.plot(Qx,Qy)
#ad.set_xlim(-numpy.pi,numpy.pi)
#ad.set_ylim(-numpy.pi,numpy.pi)

pylab.show(block=True)

quit()
