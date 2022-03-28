import pandas as pd
import matplotlib.pyplot as plt

#Load file
dt=pd.read_csv("sevn_output/output_0.csv")
#Give a look to the columns
print(dt.columns)
#Consider only the final states
dt=dt.drop_duplicates(["ID","name"], keep='last')

#Create filter indexes
idx0 = (dt.RemnantType_0==6)
idx1 = (dt.RemnantType_1==6)
idxb0 = idx0  & dt.Semimajor.notnull()
idxb1 = idx1  & dt.Semimajor.notnull()
idxm0 = idxb0 & (dt.GWtime + dt.BWorldtime  <= 14000)
idxm1 = idxb1 & (dt.GWtime + dt.BWorldtime  <= 14000)

#Filter and join masses in new dataframe
AllBH = pd.concat([dt[idx0].Mass_0,dt[idx1].Mass_1])
BoundBH = pd.concat([dt[idxb0].Mass_0,dt[idxb1].Mass_1])
MergingBH = pd.concat([dt[idxm0].Mass_0,dt[idxm1].Mass_1])


#PLot
plt.hist(AllBH,histtype="step",label="All",lw=3)
plt.hist(BoundBH,histtype="step",label="Bound",lw=3)
plt.hist(MergingBH, histtype="step",label="Merging",lw=3)
plt.yscale("log")
plt.legend(fontsize=18)
plt.xlabel("BH mass [M$_\odot$]",fontsize=18)
plt.ylabel("$N$",fontsize=18)
plt.gca().tick_params(axis='both', which='major', labelsize=18)
plt.tight_layout()
plt.savefig("example1.pdf")
plt.show()
