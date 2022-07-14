import numpy as np
import matplotlib.pyplot as plt
from math import *


def readFile(fileName):
        fileObj = open(fileName, "r") #opens the file in read mode
        mzenerg = fileObj.readline().strip().split(" ") #puts the file into an array
        info=fileObj.readline().strip();
        fileObj.close()
        return mzenerg,info
    
res,title=readFile("../data/mz_integral.txt");


res=np.array(res);
res=res.astype("float");
mzenerg=res[::2];
mztim=res[1::2];
fig1, ax1 = plt.subplots(figsize = (14, 14));
plt.plot(mztim,mzenerg,label='E_dm');
plt.title(title);
plt.savefig("../figs/"+title+".png", dpi=100)