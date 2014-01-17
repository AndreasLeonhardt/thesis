
from numpy import *
from pylab import *
from subprocess import call
import matplotlib.pyplot as plt
import array
import os

call("./StaggeredMagnetization-build-Desktop_Qt_5_0_0_GCC_64bit_SDK-Release/StaggeredMagnetization")

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
