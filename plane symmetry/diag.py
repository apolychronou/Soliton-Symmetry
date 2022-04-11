import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from math import cos,sin,pi


def cartesian(mr,mphi,mz,phi):
    m1=np.zeros([len(mr),len(phi)]);
    m2=np.zeros([len(mr),len(phi)]);
    m3=np.zeros([len(mr),len(phi)]);

    for i in range(0,len(m1)):
        for j in range(0,len(phi)):
            m1[i][j]=cos(phi[j])*mr[i]-sin(phi[j])*mphi[i];
            m2[i][j]=sin(phi[j])*mr[i]+cos(phi[j])*mphi[i];
            m3[i][j]=mz[i];
    return m1,m2,m3

def cartesian_coord(r,phi):
    x=np.zeros([len(r),len(phi)]);
    y=np.zeros([len(r),len(phi)]);
    for i in range(0,len(r)):
        for j in range(0,len(phi)):
            x[i][j]=r[i]*cos(phi[j]);
            y[i][j]=r[i]*sin(phi[j]);
    return x,y
def readFile(fileName):
        fileObj = open(fileName, "r") #opens the file in read mode
        energ=fileObj.readline();
        words = fileObj.read().split(" ") #puts the file into an array
        fileObj.close()
        return words,energ

phi=pi/2;
a,energy_info=readFile("./data/skyrm.txt");
a=a[0:-1];
j=int(len(a)/3);
r=np.zeros(j);
for i in range(0,len(r)):
    r[i]=i;
r=0.1*r;
mr=np.array(a[0::3]);
mr=mr.astype(float);
mphi=np.array(a[1::3]);
mphi=mphi.astype(float);
mz=np.array(a[2::3]);
mz=mz.astype(float);

m=np.array(a);
m=m.astype(float);

plt.plot(r,mz,'g',label="m_z");
plt.xticks(np.arange(min(r), max(r), 1));
plt.title("m vs r",fontsize=30);




plt.plot(r,mr,'b',label="m_r");
plt.plot(r,mphi,'r',label="m_phi");
plt.legend();

figure = plt.gcf()
figure.set_size_inches(14, 14);
matplotlib.rc('xtick', labelsize=20); 
matplotlib.rc('ytick', labelsize=20);

# plt.savefig("./figs/phi90.png", dpi=100);

# x=np.zeros(len(r));
# y=np.zeros(len(r));
# for i in range(0,len(r)):
#     x[i]=mr[i]*cos(mphi[i]);
#     y[i]=mr[i]*sin(mphi[i]);

# for i in range(0,len(r)):
#     print((x[i]**2+y[i]**2+mz[i]**2));

s=80;
# print("cylindrical")
phi=np.linspace(0,2*pi,len(r[:int(s/3):]))
m1,m2,m3=cartesian(mr[:s:],mphi[:s:],mz[:s:],phi);
X,Y=cartesian_coord(r[:s], phi);
l=7;
fig, ax = plt.subplots(figsize = (14, 12));
plt1=ax.quiver(X[::l], Y[::l], m1[::l], m2[::l],m3[::l],scale=4.5,headlength=4,\
          headwidth=4,cmap=plt.cm.jet,units='xy');
plt.colorbar(plt1,ax=ax)
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20)

# ax.set_aspect('equal')
# ticks=np.arange(X[::l].min(),X[::l].max(),0.5);
# ax.xaxis.set_ticks(ticks)
# ax.yaxis.set_ticks(ticks)

# plt.savefig("./figs/x-y-phi90.png", dpi=100);

# for i in range(0,len(r)):
#     print((mr[i]**2+mphi[i]**2+mz[i]**2));


# for i in range(0,len(r)):
#     print((m1[i]**2+m2[i]**2+mz[i]**2));

for i in range(0,len(m),3):
    norm=m[i]**2+m[i+1]**2+m[i+2]**2;
    # print(norm);
    if abs(1-norm)>1e-1:
        print(norm);
        
print(energy_info);
