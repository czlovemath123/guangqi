import sys
sys.path.insert(1,'../../scripts')
from assemble_1d_data import *
from math import *
import matplotlib
import matplotlib.pyplot as plt
import os.path
plt.rcParams.update({'font.size': 16})
np.set_printoptions(threshold=sys.maxsize)

def read_lc(filename):
    lc=np.zeros(402)
    f=open(filename,'r')
    for i in range(402):
        line=f.readline()
        a=[float(s) for s in line.split()]
        lc[i]=a[0]
    return lc

cwd=os.getcwd()
os.chdir(cwd+'/lc')
lc_nocol1=read_lc('lc_nocol_v1.01.data')
lc_nocol2=read_lc('lc_nocol_v1.1.data')
lc_nocol3=read_lc('lc_nocol_v1.2.data')
#lc_nocol4=read_lc('lc_nocol_v1.3.data')
#lc_nocol5=read_lc('lc_nocol_v1.4.data')
lc_col1=read_lc('lc_col_v1.01.data')
lc_col2=read_lc('lc_col_v1.1.data')
lc_col3=read_lc('lc_col_v1.2.data')
#lc_col4=read_lc('lc_col_v1.3.data')
#lc_col5=read_lc('lc_col_v1.4.data')
t0=np.linspace(0,4e6,num=403)
t=np.zeros(402)
t=t0[1:403]/day

plt.subplots(figsize=(12,6))
ln1=plt.plot(t,lc_nocol1,label='1.01vesc')
ln2=plt.plot(t,lc_nocol2,label='1.1vesc')
ln3=plt.plot(t,lc_nocol3,label='1.2vesc')
#ln4=plt.plot(t,lc_nocol4,label='1.3vesc')
#ln5=plt.plot(t,lc_nocol5,label='1.4vesc')
plt.legend()
plt.xlim(t[0],t[401])
plt.ylim(2e2,2.5e5)
plt.xlabel('day')
plt.ylabel(r'$L/L_{\odot}$')
plt.yscale('log')
plt.tight_layout()
#plt.show()
plt.savefig('nocol_lcs.png',dpi=200)
#plt.savefig('nocol_lcs.pdf',dpi=300)
plt.close()


plt.subplots(figsize=(12,6))
ln1=plt.plot(t,lc_col1,label='1.01vesc')
ln2=plt.plot(t,lc_col2,label='1.1vesc')
ln3=plt.plot(t,lc_col3,label='1.2vesc')
#ln4=plt.plot(t,lc_col4,label='1.3vesc')
#ln5=plt.plot(t,lc_col5,label='1.4vesc')
plt.legend()
plt.xlim(t[0],t[401])
plt.ylim(2e2,2.5e5)
plt.xlabel('day')
plt.ylabel(r'$L/L_{\odot}$')
plt.yscale('log')
plt.tight_layout()
#plt.show()
plt.savefig('col_lcs.png',dpi=200)
#plt.savefig('col_lcs.pdf',dpi=300)
plt.close()
