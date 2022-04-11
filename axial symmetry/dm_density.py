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
    
header,rows=readFile("./data/density.csv");

r=[];
e_dm=[];
# e_ex=[];
# e_an=[];
for i in rows:
    r.append(i[header.index("r")]);
    e_dm.append(i[header.index("E_dm")]);
    # e_ex.append(i[header.index("E_ex")]);
    # e_an.append(i[header.index("E_an")]);

e_dm=np.array(e_dm);
e_dm=e_dm.astype("float");

# e_an=np.array(e_an);
# e_an=e_an.astype("float");

# e_ex=np.array(e_ex);
# e_ex=e_ex.astype("float");

r=np.array(r);
r=r.astype("float");

fig1, ax1 = plt.subplots(figsize = (14, 14));
plt.plot(r,e_dm,label='E_dm');
# plt.plot(a,e_ex,label='E_ex');
# plt.plot(a,e_an,label='E_an');
ax1.legend();
plt.xlabel('r',fontsize=50);
plt.ylabel('DM',fontsize=50);  
matplotlib.rc('xtick', labelsize=30); 
matplotlib.rc('ytick', labelsize=40);
# plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0));
plt.xticks(np.arange(0,r.max()+4,0.5));
plt.title("DM_Density vs r",fontsize=50);
plt.xlim(xmin=0,xmax=5)
plt.ylim(ymax=10);