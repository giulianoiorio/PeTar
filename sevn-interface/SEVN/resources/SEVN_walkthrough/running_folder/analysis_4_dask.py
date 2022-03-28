import pandas as pd
import dask.dataframe as dd
import matplotlib.pyplot as plt
import numpy as np

#Load file
dt=dd.read_csv("sevn_output/output_*.csv")
#Filter BHBH binaries
#Notice we don't have to drop duplicates, since the simulation
#stop when both stars are remnant, so we are sure that with this condition
#we are taking only one row per system
idxb = (dt.RemnantType_0==6) & (dt.RemnantType_1==6) & dt.Semimajor.notnull() & (dt.BWorldtime + dt.GWtime <14000)
dt=dt[idxb]

#Load evolved file
dte=dd.read_csv("sevn_output/evolved_*.dat",sep='\s+')
dte=dte.rename(columns={'#ID': 'ID','Mass_0':"Mzams_0", 'Mass_1':"Mzams_1"})
#Join the two dataset
dt = dt.merge(dte, on=["ID","name"],  how="inner", suffixes=("","_ini") )


#CE from LOGFILE
import re
#regex string
regex_str=r'B;\d+;(\d+);CE;'
#loop over logfiles
import glob
logfiles = glob.glob("sevn_output/logfile*")
CE_array = []
for logfile in logfiles:
    with open(logfile,"r") as fo:
        ma=re.findall(regex_str,fo.read())
    #transform to a dataframe with only one column ID
    BID = pd.DataFrame({"ID":np.asarray(ma,dtype=int)})
    #Group by ID and count, this are the number of NCE occurences for each system
    BID = BID.groupby(["ID"]).size().reset_index(name="NCE")
    CE_array.append(BID)
#concatenate
CEdf = pd.concat(CE_array)
print(CEdf)
#The CEdf contains the ID of al lthe binaries that go through at least 1 CE,
#so in order to find systems that not trigger CEs we can just perform a left join
# (everything from the left, just stuff in common from the right)
#that is we join everything
dt=dt.merge(CEdf,on=['ID'],how="left")
print(dt.columns)
#then the systems no included in CEdf will have a null value in the column NCE,
#we transfotm thid to 0
dt["NCE"]=dt["NCE"].fillna(0)
#transform do pandas array
df=dt.compute()

#Plot
plt.figure(figsize=(10,5))

plt.subplot(1,2,1)
plt.scatter(df.Mzams_0,df.Mass_0,c=df.NCE,cmap="jet",vmin=0,vmax=3,edgecolor="k",s=30)
plt.scatter(df.Mzams_1,df.Mass_1,c=df.NCE,cmap="jet",vmin=0,vmax=3,edgecolor="k",s=30)
plt.xscale("log")
plt.yscale("log")
plt.ylabel("BH mass [M$_\odot$]",fontsize=18)
plt.xlabel("$M\mathrm{zams}$  [M$_\odot$]",fontsize=18)
plt.gca().tick_params(axis='both', which='major', labelsize=18)

plt.subplot(1,2,2)
plt.scatter(df.a,df.Mass_0,c=df.NCE,cmap="jet",vmin=0,vmax=3,edgecolor="k",s=30)
plt.scatter(df.a,df.Mass_1,c=df.NCE,cmap="jet",vmin=0,vmax=3,edgecolor="k",s=30)
cbar=plt.colorbar(pad=0)
cbar.ax.tick_params(axis='both', which='major', labelsize=16)
cbar.set_label(label="$N$ Common Envelope",size=15)

plt.xlabel("Semimajor initial  [R$_\odot$]",fontsize=18)
plt.ylabel("BH mass [M$_\odot$]",fontsize=18)
plt.gca().tick_params(axis='both', which='major', labelsize=18)
plt.xscale("log")
plt.yscale("log")

plt.tight_layout()
plt.savefig("analysis4.png")
plt.show()
