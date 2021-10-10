from math import *
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
np.set_printoptions(threshold=sys.maxsize)

def mdot(t):
    if (t<2):
        mdot=2.5*(1.0-0.1*t)
    else:
        mdot=0
    return mdot
def vejecta(t):
    if (t<2):
        vejecta=3*(1.0-0.1*t)
    else:
        vejecta=0
    return vejecta

t=np.linspace(0,10,101)
mdot_t=np.zeros(101)
vejecta_t=np.zeros(101)
for i in range(101):
    mdot_t[i]=mdot(t[i])
    vejecta_t[i]=vejecta(t[i])
print(mdot_t)

plt.figure(figsize=(12,8))
plt.plot(t,mdot_t,label='massloss')
plt.legend()
plt.xlim(0,10)
plt.ylabel(r'$M_{\odot}\cdot yr^{-1}$')
plt.xlabel('day')
plt.savefig('massloss.pdf',dpi=300)
#plt.show()
plt.close()

plt.figure(figsize=(12,8))
plt.plot(t,mdot_t,label='vejecta')
plt.legend()
plt.xlim(0,10)
plt.ylabel(r'$v_{esc}$')
plt.xlabel('day')
plt.savefig('vejecta.pdf',dpi=300)
#plt.show()
plt.close()
