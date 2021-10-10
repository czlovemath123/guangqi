from numpy import *
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
f=open('something.txt')
loc=[]
for i in range(400):
    line=f.readline()
    a=[float(s) for s in line.split()]
    loc.append(a)
loc=transpose(asarray(loc))
c1=plt.Circle((0,0),1.0,color='b',fill=False)
c2=plt.Circle((0,0),2.6,color='g',fill=False)
x=loc[0]
y=loc[1]
ax=plt.subplot()
ax.add_artist(c1)
ax.add_artist(c2)
plt.plot(x,y,'r.')
plt.plot(0,2,'ko')
plt.xlim(-3,3)
plt.ylim(-3,3)
ax.set_aspect(1)
plt.show()
