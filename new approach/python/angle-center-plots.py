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
    
header,rows=readFile("../data/center.csv");

a=[];
e_dm=[];
e_ex=[];
e_an=[];
for i in rows:
    a.append(i[header.index("a")]);
    e_dm.append(i[header.index("E_dm")]);
    e_ex.append(i[header.index("E_ex")]);
    e_an.append(i[header.index("E_an")]);

e_dm=np.array(e_dm);
e_dm=e_dm.astype("float");

e_an=np.array(e_an);
e_an=e_an.astype("float");

e_ex=np.array(e_ex);
e_ex=e_ex.astype("float");

a=np.array(a);
a=a.astype("float");

fig1, ax1 = plt.subplots(figsize = (14, 14));
plt.plot(a,e_dm,label='E_dm');
plt.plot(a,e_ex,label='E_ex');
plt.plot(a,e_an,label='E_an');
ax1.legend();
plt.xlabel('a',fontsize=50);
plt.ylabel('Energies',fontsize=50);  
matplotlib.rc('xtick', labelsize=20); 
matplotlib.rc('ytick', labelsize=40);
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0));
plt.xticks(np.arange(0,a.max()+4,4));
plt.title("Energies vs a",fontsize=50);
plt.savefig("./figs/energies-vs-a.png");

header,rows=readFile("./data/angle.csv");

phi=[];
e_dm=[];
e_ex=[];
e_an=[];
for i in rows:
    phi.append(i[header.index("phi")]);
    e_dm.append(i[header.index("E_dm")]);
    e_ex.append(i[header.index("E_ex")]);
    e_an.append(i[header.index("E_an")]);

e_dm=np.array(e_dm);
e_dm=e_dm.astype("float");

e_an=np.array(e_an);
e_an=e_an.astype("float");

e_ex=np.array(e_ex);
e_ex=e_ex.astype("float");

phi=np.array(phi);
phi=phi.astype("float");

fig1, ax1 = plt.subplots(figsize = (14, 14));
plt.plot(phi,e_dm,label='E_dm');
ax1.legend();
plt.xlabel('phi',fontsize=40);
plt.ylabel('e_dm',fontsize=40);  
matplotlib.rc('xtick', labelsize=20); 
matplotlib.rc('ytick', labelsize=40);
# plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0));
plt.xticks(np.arange(phi.min(),phi.max()+0.5,0.5));
plt.title("e_dm vs phi",fontsize=50);
plt.savefig("./figs/edm-vs-phi.png");

