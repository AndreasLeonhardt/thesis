import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import array
import os

plt.rcParams.update({'font.size':22})

filename = sys.argv[1]
A=loadtxt(filename, skiprows=0)

nx=int(A[-1,0])
ny=int(A[-1,0])

X=np.zeros((nx+1,ny+1))
Y=np.zeros((nx+1,ny+1))
Em=np.zeros((nx+1,ny+1))
Ep=np.zeros((nx+1,ny+1))
eps=np.zeros((nx+1,ny+1))

for i in range(nx+1):
	for j in range(ny+1):
		X[i,j]=(A[i*(ny+1),2])/(2*np.pi)
		Y[i,j]=(A[j,3])/(2*np.pi)
		Em[i,j]=(A[i+j*(nx+1),5])
		Ep[i,j]=(A[i+j*(nx+1),6])
		eps[i,j]=A[i+j*(nx+1),4]#*.259#eV

for j in range(ny):
	X[nx,j]=-X[0,j]
	Ep[nx,j]=Ep[0,j]
	Em[nx,j]=Em[0,j]
	eps[nx,j]=eps[0,j]
	Y[nx,j]=Y[0,j]

for i in range(nx+1):
	X[i,ny]=X[i,0];
	Y[i,ny]=-Y[i,0];
	eps[i,ny]=eps[i,0];
	Em[i,ny]=Em[i,0]
	Ep[i,ny]=Ep[i,0]

print(-np.max(Ep)+np.min(Em))

# reduced zone
nxh=(nx+1)/2
nyh=(ny+1)/2

rX = np.zeros((nxh,nyh))
rY = np.zeros((nxh,nyh))
rEp= np.zeros((nxh,nyh))
rEm= np.zeros((nxh,nyh))

for i in range(nxh):
	for j in range(nyh):
		rX[i,j]=X[nxh-i+j,i+j]
		rY[i,j]=Y[nxh-i+j,i+j]
		rEp[i,j]=Ep[nxh-i+j,i+j]
		rEm[i,j]=Em[nxh-i+j,i+j]

pTmax= 7*nxh/2
ix=3*nxh/2-1
iy=3*nyh/2-1
it=0
pT=np.zeros(pTmax)
pEp=np.zeros(pTmax)
pEm=np.zeros(pTmax)

for i in range(pTmax):
	pT[i] = it
	pEp[i]=Ep[ix,iy]
	pEm[i]=Em[ix,iy]
	if i<nxh/2:
		ix+=+1
		iy+=-1
		it+=np.sqrt(2)
	elif i<3*nxh/2:
		iy+=1
		it+=1
	elif i<5*nxh/2:
		ix+=-1
		iy+=-1
		it+=np.sqrt(2)
	else: 
		ix+=-1
		it+=1
		 


t=np.arange(np.floor(eps.min()),np.ceil(eps.max()),1.)

f = figure()

#ad = f.add_subplot(111,projection='3d')
#ad.plot_wireframe(rX,rY,rEm )
#ad.plot_wireframe(rX,rY,rEp,color='red')
#xlabel('$k_x$')
#ylabel('$k_y$')
#zlabel('E/t')


#ad = f.add_subplot(111)
#ad.arrow(-.405,-.095,.5-0.04/np.sqrt(2),.5-0.04/np.sqrt(2),head_width=.02,head_length=.04,fc="k", ec="k" )
#p=ad.contour(X,Y,eps,t)
#cb=f.colorbar(p,ax=ad)
#plt.xticks(np.arange(-.5,.6,.25))
#plt.yticks(np.arange(-.5,.6,.25))
#xlabel('j')
#ylabel('E')

ad = f.add_subplot(111)
ad.plot(pT,pEp,color='red')
ad.plot(pT,pEm,color='blue')
xlabel('$\\vec{q}$')
ylabel('$\\frac{E}{t}$')
xticks([0,nxh/2*np.sqrt(2),(2+np.sqrt(2))*nxh/2,(2+2*np.sqrt(2))*nxh/2,(2+3*np.sqrt(2))*nxh/2,(4+3*np.sqrt(2))*nxh/2],['S','X','M','S','$\Gamma$','X'])
ylim(-6,5)

show(block=True)

quit()
