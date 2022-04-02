import numpy as np
import matplotlib.pyplot as plt
from math import *

FILES=1;
START=0
kinisi=np.zeros([FILES,2],type(int));
for li in range(START,FILES):
    def readFile(fileName):
            fileObj = open(fileName, "r") #opens the file in read mode
            fileObj.readline();
            words = fileObj.read().split(" ") #puts the file into an array
            fileObj.close()
            return words
    
    
    # m=readFile("./data/vort"+str(li)+".txt");
    m=readFile("./data/3D"+str(li)+".txt");
    m=m[0:-1];
    N=200;
    NZ=100;
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
    mz_row=mz;
    mz=mz.reshape(NZ,N);
    mr=mr.reshape(NZ,N);
    mphi=mphi.reshape(NZ,N);
    m=np.array(m);
    m=m.astype(float);
    m=m.reshape(NZ,N,3);
    

# ------------------- CONTOUR PLOT ------------------
    
    # fig = plt.figure()

    # R, Z = np.meshgrid(r, z);
    # plt.contour(R, Z, mz,levels=[-0.8,-0.4,0,0.4,0.8], colors='black');


    # plt.xlim([0,15]);
    # plt.ylim([-10,10]);
    # plt.savefig("./figs/contour_t="+str(li)+".png", dpi=100)

#  ----------------- 3D PLOT ----------------------

    R, Z = np.meshgrid(r, z);
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_wireframe(R, Z, mz);
    ax.set_xlabel('r')
    ax.set_ylabel('z')
    ax.set_zlabel('m_z');
    
    
    
    fig.set_size_inches(14, 12)
    plt.savefig("./figs/mz"+str(li)+".png", dpi=100)
    
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_wireframe(R, Z, mr);
    ax.set_xlabel('r')
    ax.set_ylabel('z')
    ax.set_zlabel('m_r');
    fig.set_size_inches(14, 12)
    plt.savefig("./figs/mr"+str(li)+".png", dpi=100)
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_wireframe(R, Z, mphi);
    ax.set_xlabel('r')
    ax.set_ylabel('z')
    ax.set_zlabel('m_phi');
    fig.set_size_inches(14, 12);
    plt.savefig("./figs/mphi"+str(li)+".png", dpi=100)

    
    result_vort8 = np.where(m == mz.min());
    
    indices = np.where(result_vort8[2]==2)[0][0]
    lista=[];

    for i in result_vort8:
        lista.append(i[indices]);
        
    # lilist.append(lista); 
    kinisi[li]=np.array(lista[0:2]);
kinisi=kinisi.astype(int);

out=np.zeros([FILES,2]);
j=0;
for i in kinisi:
    out[j][0]=z[i[0]]
    out[j][1]=r[i[1]];
    j=j+1;
np.savetxt("movement.csv", out, delimiter=",");
    # plt.savefig("./figs/mphi.png", dpi=100)

# fig = plt.figure()
# plt.plot(r[:],mz[0][:],'g',label="m_z");
# plt.xticks(np.arange(r[0],r[-1],0.5))
# plt.title("m vs r");




# plt.plot(r[:],mr[0][:],'b',label="m_r");
# plt.plot(r[:],mphi[0][:],'r',label="m_phi");
# plt.legend();

# figure = plt.gcf()
# figure.set_size_inches(14, 12)
# plt.savefig("./figs/skyrmion.png")


# for i in range(0,NZ):
#     if(np.linalg.norm(mr[i][:]-mr[0][:])>1e-1):
#         print("diff mr\n");
#         print(np.linalg.norm(mr[i][:]-mr[0][:]));
#     if(np.linalg.norm(mphi[i][:]-mphi[0][:])>1e-1):
#         print("diff mphi\n");
#         print(np.linalg.norm(mphi[i][:]-mphi[0][:]));
#     if(np.linalg.norm(mz[i][:]-mz[0][:])>1e-1):
#         print("diff mz\n");
#         print(np.linalg.norm(mz[i][:]-mz[0][:]));


# for i in range(0,N):
#     for j in range(0,NZ):
#         norm=mr[j][i]**2+mz[j][i]**2+mphi[j][i]**2;
#         if abs(norm-1)>1e-8:
#             print(norm);
