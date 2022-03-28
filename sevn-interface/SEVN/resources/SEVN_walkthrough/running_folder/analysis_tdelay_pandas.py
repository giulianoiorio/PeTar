import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#Load file
dt=pd.read_csv("sevn_output/output_0.csv")
#Give a look to the columns
print(dt.columns)
#Consider only the final states
dt=dt.drop_duplicates(["ID","name"], keep='last')

#Create filter indexes
idx_bound = (dt.RemnantType_0==6) &  (dt.RemnantType_1==6) & dt.Semimajor.notnull()
dt =dt[idx_bound]

#Create new columnd
dt["tdelay"] = dt.BWorldtime + dt.GWtime



#PLot
bins=np.logspace(np.log10(100),np.log10(1.4e4),50)
plt.hist(dt.tdelay,range=(1e-2,1e5),histtype="step",label="All",lw=3,bins=bins)
#plt.hist(BoundBH,histtype="step",label="Bound",lw=3)
#plt.hist(MergingBH, histtype="step",label="Merging",lw=3)
#plt.yscale("log")
#plt.legend(fontsize=18)
#plt.xlabel("BH mass [M$_\odot$]",fontsize=18)
#plt.ylabel("$N$",fontsize=18)
#plt.gca().tick_params(axis='both', which='major', labelsize=18)
#plt.tight_layout()
#plt.savefig("example1.pdf")
plt.xscale("log")
plt.show()
