import numpy as np
import matplotlib.pyplot as plt
from math import *

FILES=11;
START=0
kinisi=np.zeros([FILES,2],type(int));
minim=np.zeros(FILES);

N=300;
NZ=200;
dr=0.1;
dz=0.1;
centz=NZ/2;
centr=40;
energy_info=[]
for li in range(START,FILES):
    def readFile(fileName):
            fileObj = open(fileName, "r") #opens the file in read mode
            words = fileObj.readline().split(" ") #puts the file into an array
            energy_info=fileObj.readline();
            fileObj.close()
            return words,energy_info
    
    
    # m=readFile("./data/vort"+str(li)+".txt");
    m,e=readFile("./data/vort"+str(li)+".txt");
    m=m[0:-1];
    energy_info.append(e);
    r=np.zeros(N);
    z=np.zeros(NZ);
    for i in range(0,len(r)):
        r[i]=i;
    for i in range(0,len(z)):
        z[i]=-NZ/2*dz+i*dz;
    r=r*dr;
    
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
    fig = plt.figure()
    fig.set_size_inches(14, 14)

    R, Z = np.meshgrid(r, z);
    plt.contour(R, Z, mz,levels=[-0.9,-0.5,0,0.9], colors='black');


    plt.xlim([0,10]);
    plt.ylim([-7,7]);
    plt.title("t="+str(li));
    # plt.savefig("./figs/contour_t="+str(li)+".png", dpi=100)

#  ----------------- 3D PLOT ----------------------

    # R, Z = np.meshgrid(r, z);
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.plot_wireframe(R, Z, mz);
    # ax.set_xlabel('r')
    # ax.set_ylabel('z')
    # ax.set_zlabel('m_z');
    
    
    
    # fig.set_size_inches(14, 12)
    # plt.savefig("./figs/mz"+str(li)+".png", dpi=100)
    
    
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.plot_wireframe(R, Z, mr);
    # ax.set_xlabel('r')
    # ax.set_ylabel('z')
    # ax.set_zlabel('m_r');
    # fig.set_size_inches(14, 12)
    # plt.savefig("./figs/mr"+str(li)+".png", dpi=100)
    
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.plot_wireframe(R, Z, mphi);
    # ax.set_xlabel('r')
    # ax.set_ylabel('z')
    # ax.set_zlabel('m_phi');
    # fig.set_size_inches(14, 12);
    # plt.savefig("./figs/mphi"+str(li)+".png", dpi=100)


    #------- movement coordinates ----------- 
    plz1=centz-15;
    plz2=centz+15;
    plr1=centr-5;
    plr2=centr+5;
    while(plz1<0):
        plz1+=1;
        plz2+=1;
    while(plr1<0):
        plr1+=1;
        plr2+=1;
    plz1=int(plz1);plz2=int(plz2);
    plr1=int(plr1);plr2=int(plr2);
  
    result_vort8 = np.where(mz== mz[plz1:plz2,plr1:plr2].min());
    
    # indices = np.where(result_vort8[2]==2)[0][0]
    # indices = result_vort8[0][0]
    lista=[];
    minim[li]=mz[plz1:plz2,plr1:plr2].min();

    cr1=0;
    cr2=0;
    for i in result_vort8:
        # lista.append(i[indices]);
        cr1+=1;
        if cr1==1:
            for k in i:
                if(k>=plz1 and k<=plz2):
                    lista.append(k);
                else:
                    cr2+=1;
        else:
            lista.append(i[cr2]);
            


        
    # lilist.append(lista); 
    kinisi[li]=np.array(lista[0:2]);
    centz=int(kinisi[li][0]);
    centr=int(kinisi[li][1]);
    
    
kinisi=kinisi.astype(int);

out=np.zeros([FILES,3]);
j=0;
for i in kinisi:
    out[j][0]=z[i[0]]
    out[j][1]=r[i[1]];
    out[j][2]=minim[j];
    j=j+1;
np.savetxt("./figs/movement.csv", out, delimiter=",");
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
