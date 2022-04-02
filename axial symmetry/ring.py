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
        fileObj.readline();
        words = fileObj.read().split(" ") #puts the file into an array
        fileObj.close()
        return words

m=readFile("./data/vortex.txt");
m=m[0:-1];
N=600;
NZ=300;
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
# mz_row=mz;
mz=mz.reshape(NZ,N);
mr=mr.reshape(NZ,N);
mphi=mphi.reshape(NZ,N);
m=np.array(m);
m=m.astype(float);
# mr[np.where(abs(mr)<1e-3)]=0;
# mphi[np.where(abs(mr)<1e-3)]=0;
# mz[np.where(abs(mr)<1e-3)]=0;



# m=m.reshape(NZ,N,3);

s1=0;
s2=70;
phi=np.linspace(0,2*pi,len(r[s1:int(s2/2):]));
X,Y=cartesian_coord(r[s1:s2], phi);
i=0;
fig, ax = plt.subplots(figsize = (16, 14));
fig2, ax2 = plt.subplots(figsize = (16, 14));

for z_slice in range (145,155,1):
    i=i+1;
    
    # --------- 2D PLOTS ---------------------
    
    m1,m2,m3=cartesian(mr[z_slice][s1:s2:],mphi[z_slice][s1:s2:],mz[z_slice][s1:s2:],phi);
    
    l=7;
    
    
    fig, ax = plt.subplots(figsize = (16, 14));
    plt1=ax.quiver(X[::l], Y[::l], m1[::l], m2[::l],m3[::l],scale=4,headlength=4,\
              headwidth=4,cmap=plt.cm.jet,scale_units='xy');
    ticks=np.arange(-r[s2],r[s2]+0.5,0.5);
    ax.xaxis.set_ticks(ticks)
    ax.yaxis.set_ticks(ticks)
    plt.colorbar(plt1,ax=ax, cmap=plt.cm.jet)

    plt.title("z_slice="+str(round(z[z_slice],2)));
    
    
    # plt.savefig("./figs/z_slice="+str(round(z[z_slice],2))+".png", dpi=100)
    
   # ---------- 1D PLOTS ----------------------
    
    # -------- mr-mphi
#     l=5;
#     y=np.zeros(len(r[s1:s2:l]));
#     y=y-z[z_slice];
#     # y=np.ones(len(r[s1:s2:l]));
#     # fig, ax = plt.subplots(figsize = (18, 10));
#     plt1=ax.quiver(r[s1:s2:l],y,mr[z_slice][s1:s2:l],mphi[z_slice][s1:s2:l],\
#               mz[z_slice][s1:s2:l],cmap=plt.cm.jet,scale=25,headlength=4,headwidth=4);
#     # plt.title("z_slice="+str(round(z[z_slice],2)));
#     title = ax.set_title("mr-mphi", fontsize='large')

    
#     # plt.savefig("./figs/z_slice="+str(round(z[z_slice],2))+".png", dpi=100)
#     # plt.savefig("./figs/"+str(i)+".png", dpi=100)

#     # ------- mr- mz
#     plt2=ax2.quiver(r[s1:s2:l],y,mr[z_slice][s1:s2:l],mz[z_slice][s1:s2:l],\
#               mphi[z_slice][s1:s2:l],cmap=plt.cm.jet,scale=25,headlength=4,headwidth=4);

#     title = ax2.set_title("mr-mz", fontsize='large')
#     # plt.title("z_slice1="+str(round(z[z_slice],2)));
    
#     #plt2.savefig("./figs/z_slice="+str(round(z[z_slice],2))+".png", dpi=100)
#     #plt2.savefig("./figs/"+str(i)+".png", dpi=100)

# plt.colorbar(plt1,ax=ax, cmap=plt.cm.jet)
# plt.colorbar(plt2,ax=ax2, cmap=plt.cm.jet)



