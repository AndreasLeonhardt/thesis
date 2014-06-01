# python script for plotting results from Susceptibility
# call with data file as argument
import  numpy as np
import pylab
from subprocess import call
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import array
import os

plt.rcParams.update({'font.size':26})

# load data from file, using the first argument as filename.
filename = pylab.sys.argv[1]
A=pylab.loadtxt(filename, skiprows=2)

# get max coordinates (q,w)
# (might need modification for offset, e.g. not starting with w=0.)
tmax=int(A[-1,0])+1
n_w =int(len(A[:,0])/(tmax))

# specifier for type of plot, as written from Susceptibility.
# where b indicates a bar (q+Q) and r and i the real and imaginary parts respectiviely

qx=1
qy=2
n_xr=4
n_xi=5
n_xbr=6
n_xbi=7
n_yr=8
n_yi=9
#n_ybr=10
#n_ybi=11
n_z1r=10
n_z1i=11
#n_z1br=14
#n_z1bi=15
n_z2r=12
n_z2i=13
#n_z2br=18
#n_z2bi=19


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
N_x = abs(np.sqrt(N))
lmb=U/N
m=float(params[7])

# create np arrays for coordinates and values
Val=np.ndarray((tmax,n_w))          
T=np.array(range(tmax)).astype(float)
W=np.array(range(n_w )).astype(float)
Qx=[]
Qy=[]
Heis=[]
Heisttt=[]
t1=1.0
t2=0.2326
t3=0.1163

theta = .191986 # =11 degree, angle of octahedral rotations (e.g. they are rotated \pm\theta, their relative angle is  2*theta)



# fill with corresponding coordinates/values
pos= 0.0
for t in range (0,tmax):
	#gradually rescale T in diagonal steps (divide by sqrt(2) and shift sum accordingly)
	# diagonal in (0,1); (3,5) in terms of N_x/4=N/8
	# straight in (1,3); (5,7)
	T[t]=pos  #A[t*n_w,0]/np.sqrt(2)
	if (t<N_x/4):
		pos+=np.sqrt(2)
	elif (t<3*N_x/4.):
		pos+=1.
	elif (t<5*N_x/4.):
		pos+=np.sqrt(2)
	else:
		pos+=1.;


# plotted values are selected (x,xb,y,z1,z2) /calculated (chi, denominator of chi) here.
	for w in range(0,n_w):
		tw=t*n_w+w
		x=complex(A[tw,n_xr],A[tw,n_xi])
		xb=complex(A[tw,n_xbr],A[tw,n_xbi])
		y=complex(A[tw,n_yr],A[tw,n_yi])
		#yb=complex(A[tw,n_ybr],A[tw,n_ybi])
		z1=complex(A[tw,n_z1r],A[tw,n_z1i])
		#z1b=complex(A[tw,n_z1br],A[tw,n_z1bi])
		z2=complex(A[tw,n_z2r],A[tw,n_z2i])
		#z2b=complex(A[tw,n_z2br],A[tw,n_z2bi])

		
		# calculate chi (imaginary part of the susceptibility
		
		# transversal (+-)
		Val[t,w]=((-(x+y)*(1.0+lmb*(xb+y))+lmb*(z1+z2)*(z1+z2))/((1.0+lmb*(x+y))*(1.0+lmb*(xb+y))-lmb**2*(z1+z2)*(z1+z2))*(np.cos(theta))**2+(-(xb+y)*(1.0+lmb*(x+y))+lmb*(z1+z2)*(z1+z2))/((1.0+lmb*(xb+y))*(1.0+lmb*(x+y))-lmb**2*(z1+z2)*(z1+z2))*(np.sin(theta)**2) ).imag*4
		
		# longotudinal (zz)
		temp  = 1 - lmb**2*(x-y)**2 + lmb**2*(z1-z2)**2
		tempb = 1 - lmb**2*(xb-y)**2 + lmb**2*(z1-z2)**2
		Val[t,w]+= -(( tempb*(x-y)*(1+lmb*(x-y)) - temp*lmb*(z1-z2)**2 - lmb**2*(x-xb)*(z1-z2)**2*(1+2*lmb*(x-y)))/ (temp*tempb + lmb**4*(z1-z2)**2*(x-xb)**2) ).imag*4
			
	Qx.append(A[t*n_w,1]/(2*np.pi))
	Qy.append(A[t*n_w,2]/(2*np.pi))
	Heis.append(4.0/(float(U*1))*np.sqrt(4.0-(np.cos(A[t*n_w,1])+np.cos(A[t*n_w,2]))**2)*0.258)
	
#	E=0.0
# for calculation of large U limit. Needs better expansion of \chi, taking into account that in general \varepsilon(p+Q) \ne -\varepsilon(p) 
#	for kk in range(50):
#		for ll in range(50):
#			ep=-2*t1*(np.cos(np.pi*(kk/50-1)) + np.cos(np.pi*(ll/50-1)))-4*t2*np.cos(np.pi*(kk/50-1))*np.cos(np.pi*(ll/50-1))-2*t3*(np.cos(np.pi*(kk/25-2)) + np.cos(np.pi*(ll/25-3)))
#			epq=-2*t1*(np.cos(np.pi*(kk/50-1)+A[t*n_w,1]) + np.cos(np.pi*(ll/50-1)+A[t*n_w,2]))
#			epq+= -4*t2*np.cos(np.pi*(kk/50-1)+A[t*n_w,1])*np.cos(np.pi*(ll/50-1)+A[t*n_w,2])
#			epq+= -2*t3*(np.cos(np.pi*(kk/25-2)+2*A[t*n_w,1]) + np.cos(np.pi*(ll/25-2)+2*A[t*n_w,2]))
#			E+=ep**4+epq**4 -2*ep**2*epq**2
#	E/=50*50		
#	print E
#	Heisttt.append(np.sqrt(E)/U)

# 4/(float(U))*(t-tt-ttt)*np.sqrt(4*(t-tt-ttt)**2-(t*(np.cos(A[t*n_w,1])+np.cos(A[t*n_w,2]))-2*tt*np.cos(A[t*n_w,1])*np.cos(A[t*n_w,2])-ttt*(np.cos(A[t*n_w,1]*2)+np.cos(A[t*n_w,2]*2)))**2))


for w in range(n_w):
	# scale frequency to Heisenberg model predictions, that is with 4t^2/U 
	#(since w is in units of t or t=1, that would be divided by U/t, which is U here)
	W[w]=1.18/2.0*float(A[w,3])*0.258#eV    , HB: /4*U
    
# create 2D arrays for plot function using meshgrid (see doc)
W,T=np.meshgrid(W,T)


# create figure, subfigure and plot finally
fig  = plt.figure()
#ad = fig.add_subplot(1,1,1,projection='3d')
ad= fig.add_subplot(1,1,1)

#p= ad.plot_wireframe(T,W,Val)
p=ad.pcolormesh(T,W,Val, cmap='Blues' )
#p=ad.imshow(Val.transpose(),extent = [T.min(),T.max(),W.min(),W.max()], aspect='auto',origin = 'lower',cmap = 'Blues',interpolation ='none') # RdBu, Blues
p.set_clim(-000,50000)
#p= ad.contour(T,W,Val)
#cb = fig.colorbar(p,ax=ad)

#p=ad.plot(T,Heis)

N=256/4
N=7/4*N
pylab.xlabel('$\\vec q$')
pylab.xlim([0,T.max()])
pylab.ylabel('$\\omega$ [eV]')
pylab.ylim([0,1.000])
pylab.xticks([0,np.sqrt(2)*N,(2+np.sqrt(2))*N,(2+2*np.sqrt(2))*N,(2+3*np.sqrt(2))*N,(4+3*np.sqrt(2))*N],['S','X','M','S','$\\Gamma$','X'])
#p=ad.plot(W[0,:],Val[0,:],W[96,:],Val[96,:])

# limit z axis, since results may diverge
#ad.set_zlim3d(-5000,5000)

# plot path through Brillouin zone
#ad= fig.add_subplot(1,1,1)
#p=ad.plot(Qx,Qy)
#ad.set_xlim(-.502,.502)
#ad.set_ylim(-.502,.502)
#ad.set_xlabel('$k_x/2\pi$')
#ad.set_ylabel('$k_y/2\pi$')
#p=ad.plot([.5,0,-.5,0,.5],[0,.5,0,-.5,0],':',color='0.50')
#ad.text(0,0,'$\Gamma$',horizontalalignment='right',verticalalignment='top')
#ad.text(.52,.5,'M',horizontalalignment='left',verticalalignment='top')
#ad.text(.52,0,'X',horizontalalignment='left',verticalalignment='top')
#ad.text(.25,.25,'S',horizontalalignment='right',verticalalignment='bottom')
pylab.show(block=True)

quit()
