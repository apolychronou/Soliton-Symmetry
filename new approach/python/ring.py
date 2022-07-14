import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from math import *
import matplotlib.cm as cm
from matplotlib.colors import Normalize

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

def cartesian_1D(mr,mphi):
    phi=0;
    x=np.zeros(len(mr));
    y=np.zeros(len(mr));
    for i in range(0,len(mr)):
        x[i]=cos(phi)*mr[i]-sin(phi)*mphi[i];
        y[i]=sin(phi)*mr[i]+cos(phi)*mphi[i];
    return x,y


def readFile(fileName):
        fileObj = open(fileName, "r") #opens the file in read mode
        # fileObj.readline();
        words = fileObj.readline().split(" ") #puts the file into an array
        energy=fileObj.readline();
        print(energy);
        fileObj.close()
        return words

m=readFile("../data/vort0.txt");
m=m[0:-1];
N=300;
NZ=200;
r=np.zeros(N);
z=np.zeros(NZ);
for i in range(0,len(r)):
    r[i]=i;
for i in range(0,len(z)):
    z[i]=-NZ/2*0.1+i*0.1;
r=0.1*r;

mr=np.array(m[0::3]);
mr=mr.astype(float);
mphi=np.array(m[1::3]);
mphi=mphi.astype(float);
mz=np.array(m[2::3]);
mz=mz.astype(float);
mz=mz.reshape(NZ,N);
mr=mr.reshape(NZ,N);
mphi=mphi.reshape(NZ,N);
m=np.array(m);
m=m.astype(float);


s1=0;
s2=60;
phi=np.linspace(0,2*pi,len(r[s1:int(s2/3):1]));
X,Y=cartesian_coord(r[s1:s2], phi);
i=0;
fig, ax = plt.subplots(figsize = (16, 14));
fig2, ax2 = plt.subplots(figsize = (16, 14));

start=70;
end=131;
for z_slice in range (start,end,5):
    i=i+1;
    
    # --------- 2D PLOTS ---------------------
    
    m1,m2,m3=cartesian(mr[z_slice][s1:s2:],mphi[z_slice][s1:s2:],mz[z_slice][s1:s2:],phi);
    
    l=5;
    
    
    fig1, ax1 = plt.subplots(figsize = (16, 14));
    plt3=ax1.quiver(X[::l], Y[::l], m1[::l], m2[::l],m3[::l],scale=30,headlength=3.2,\
              headwidth=3.2,cmap=plt.cm.jet,width=0.003);
    ticks=np.arange(-r[s2],r[s2]+0.5,0.5);
    ax.xaxis.set_ticks(ticks)
    ax.yaxis.set_ticks(ticks)
    plt.colorbar(plt3,ax=ax1)

    plt.title("z_slice="+str(round(z[z_slice],2)));
    
        
   # ---------- 1D PLOTS ----------------------
    
    # ##-------- mr-mphi
    l=2;
    y=np.zeros(len(r[s1:s2:l]));
    y=y+z[z_slice];
    plt1=ax.quiver(r[s1:s2:l],y,mphi[z_slice][s1:s2:l],mz[z_slice][s1:s2:l],\
              scale=25,width=0.005,pivot='middle',headlength=3.5,headwidth=3.5);
    title = ax.set_title("mr-mz", fontsize='large')
    matplotlib.rc('xtick', labelsize=30) 
    matplotlib.rc('ytick', labelsize=30)
    title = ax.set_title("mphi-mz", fontsize='large')
    ax.set_yticks(np.arange(z[start],z[end],0.2));
    ax.set_xticks(np.arange(0,r[s2]+0.5,0.5));


    

    # ##------- mr- mz---------------
    plt2=ax2.quiver(r[s1:s2:l],y,mr[z_slice][s1:s2:l],mz[z_slice][s1:s2:l],\
              scale=25,width=0.005,pivot='middle',headlength=3.5,headwidth=3.5);
    title = ax2.set_title("mr-mz", fontsize='large')
    matplotlib.rc('xtick', labelsize=30) 
    matplotlib.rc('ytick', labelsize=30)
    title = ax2.set_title("mr-mz", fontsize='large')
    ax2.set_yticks(np.arange(z[start],z[end],0.2));
    ax2.set_xticks(np.arange(0,r[s2]+0.5,0.5));






# ## -----------1 z-slice mr,mphi,mz
# fig3, ax3 = plt.subplots(figsize = (14, 14));
# plt.plot(r,mz[int(NZ/2),:],'g',label="m_z");
# plt.xticks(np.arange(min(r), max(r)+1, 1),fontsize=20);
# plt.title("m vs r",fontsize=30);
# plt.plot(r,mr[int(NZ/2),:],'b',label="m_r");
# plt.plot(r,mphi[int(NZ/2),:],'r',label="m_phi");
# plt.legend();





for i in range(0,len(m),3):
    norm=m[i]**2+m[i+1]**2+m[i+2]**2;
    # print(norm);
    if abs(1-norm)>1e-1:
        print(norm);


