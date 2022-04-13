import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib



def readFile(fileName):
        fileObj = open(fileName, "r") #opens the file in read mode
        csvreader=csv.reader(fileObj);
        header=[];
        header = next(csvreader);
        rows = []
        for row in csvreader:
            rows.append(row)
        fileObj.close();
        return header,rows
    
header,rows=readFile("../data/density.csv");

r=[];
e_dm=[];

for i in rows:
    r.append(i[header.index("r")]);
    e_dm.append(i[header.index("E_dm")]);


e_dm=np.array(e_dm);
e_dm=e_dm.astype("float");



r=np.array(r);
r=r.astype("float");

fig1, ax1 = plt.subplots(figsize = (14, 14));
plt.plot(r,e_dm,label='E_dm');
ax1.legend();
plt.xlabel('r',fontsize=50);
plt.ylabel('DM',fontsize=50);  
matplotlib.rc('xtick', labelsize=30); 
matplotlib.rc('ytick', labelsize=40);
plt.xticks(np.arange(0,r.max()+4,0.5));
plt.title("DM_Density vs r",fontsize=50);
plt.xlim(xmin=0,xmax=5)
plt.ylim(ymax=10);

plt.savefig("./figs/dm_density-vs-r.png");



## -------------- Z DENSITY

header,rows=readFile("./data/z_density.csv");
z=[];
e_dm=[];

for i in rows:
    z.append(i[header.index("z")]);
    e_dm.append(i[header.index("E_dm")]);


e_dm=np.array(e_dm);
e_dm=e_dm.astype("float");



z=np.array(z);
z=z.astype("float");

fig1, ax1 = plt.subplots(figsize = (14, 14));
plt.plot(z,e_dm,label='E_dm');
ax1.legend();
plt.xlabel('z',fontsize=50);
plt.ylabel('DM',fontsize=50);  
matplotlib.rc('xtick', labelsize=30); 
matplotlib.rc('ytick', labelsize=40);
# plt.xticks(np.arange(0,10+0.5,0.5));
plt.title("E_dm vs z",fontsize=50);
plt.xlim(xmin=0,xmax=10)
# plt.ylim(ymax=30);

plt.savefig("./figs/e_dm-vs-z.png");

