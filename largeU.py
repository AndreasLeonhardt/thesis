from sympy import *
import numpy as np


n=100
cpx,cpy,spx,spy,cqx,cqy,sqx,sqy = symbols('cpx cpy spx spy cqx cqy sqx sqy')

cp2x,cp2y,sp2x,sp2y,cq2x,cq2y,sq2x,sq2y= symbols('cp2x cp2y sp2x sp2y cq2x cq2y sq2x sq2y')

t,tt,ttt=symbols('t tt ttt')


#ep = -2*t*(cpx +cpy)
#epq = ep.subs({cpx:cpx*cqx-spx*sqx,cpy:cpy*cqy-spy*sqy})
#E=ep**4-ep**2*epq**2

#disp=0
#for i in range(n):
#	px = 2*np.pi*(i/n-.5)
#	for j in range(n):
#		py = 2*np.pi*(j/n-.5)	
#		disp+= E.subs({cpx:np.cos(px),spx:np.sin(px),cpy:np.cos(py),spy:np.sin(py)})
#disp/=n**2



ep = -2*t*(cpx +cpy)-4*tt*(cpx*cpy)-2*ttt*(cp2x+cp2y)
epq = ep.subs({cpx:cpx*cqx-spx*sqx,cpy:cpy*cqy-spy*sqy,cp2x:cp2x*cq2x-sp2x*sq2x,cp2y:cp2y*cq2y-sp2y*sq2y})


E=(ep**4-epq**4-2*ep**2*epq**2)


disp=0
for i in range(n):
	px = 2*np.pi*(i/n-.5)
	for j in range(n):
		py = 2*np.pi*(j/n-.5)	
		disp+= E.subs({cpx:np.cos(px),spx:np.sin(px),cpy:np.cos(py),spy:np.sin(py),cp2x:np.cos(2*px),sp2x:np.sin(2*px),cp2y:np.cos(2*py),sp2y:np.sin(2*py)})
disp/=n**2

print disp

