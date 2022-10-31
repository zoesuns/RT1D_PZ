import numpy as np
import matplotlib.pyplot as plt
from pyfrit.RT1D.xT_adaptive import *
import subprocess

tauTol=0.1
try:
    subprocess.call(["mkdir","2kpc"])
except:
    pass


T0=1e4
ttest=np.logspace(0,8,25)*yr2s
alphaA=(recomb_rate(T0,"HI"))
Ndot=1e57
nH=1e-2
tRec=1/(alphaA*nH)
print(tRec)
rS=(3*Ndot/(4*np.pi*alphaA*nH**2))**0.333
rI=lambda t:rS*(1-np.exp(-t/tRec))**0.333
print(rI(ttest)/Mpc2cm)
plt.plot(ttest/yr2s/1e6,rI(ttest)/Mpc2cm)
plt.xlabel(r"$\rm t[Myr]$",fontsize=20)
plt.ylabel(r"$\rm I-front position [Mpc]$")
plt.loglog()
#plt.savefig("analytical.pdf")
