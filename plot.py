
from numpy import *
from pylab import *
from subprocess import call
import matplotlib.pyplot as plt
import array
import os


filename = sys.argv[1]
A=loadtxt(filename, skiprows=0)



f = figure()
ad = f.add_subplot(111)

plot(A[:,0],A[:,1])

xlabel('U/t')
ylabel('m/2')
#title('blocking results')




show(block=True)

quit()
