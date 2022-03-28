import dask.dataframe as dd
import matplotlib.pyplot as plt
import numpy as np

#Load file
dt=dd.read_csv("sevn_output/output_*.csv")
#Give a look to the columns
print(dt.columns)
#Consider only the final states
dt=dt.drop_duplicates(["ID","name"], keep='last')

#Load evolved file
dte=dd.read_csv("sevn_output/evolved_*.dat",sep='\s+')
#Give a look to the columns
print(dte.columns)
dte=dte.rename(columns={'#ID': 'ID','Mass_0':"Mzams_0", 'Mass_1':"Mzams_1"})
#After change
print(dte.columns)

#Join the two dataset
dt = dt.merge(dte, on=["ID","name"],  how="inner", suffixes=("","_ini") )
# - on: column(s, can be a list of columns) to match during the merge of the two tables. The colum(s) has(have) to be present in both the tables
# - how: type of join to use, see documentation here and the next slide
# - suffixes: columns with the same name in the two tables (not used in on) will be renamed adding these suffixes.
#Give a look to the columns
print(dt.columns)


#Create filter indexes
idx0 = (dt.RemnantType_0==6)
idx1 = (dt.RemnantType_1==6)
idxb0 = idx0  & dt.Semimajor.notnull()
idxb1 = idx1  & dt.Semimajor.notnull()
idxm0 = idxb0 & (dt.GWtime + dt.BWorldtime  <= 14000)
idxm1 = idxb1 & (dt.GWtime + dt.BWorldtime  <= 14000)


#Filter and join masses
AllBH = dd.concat([dt[idx0].Mass_0,dt[idx1].Mass_1])
BoundBH = dd.concat([dt[idxb0].Mass_0,dt[idxb1].Mass_1])
MergingBH = dd.concat([dt[idxm0].Mass_0,dt[idxm1].Mass_1])

#Filter and join initial masses
AllBHzams = dd.concat([dt[idx0].Mzams_0,dt[idx1].Mzams_1])
BoundBHzams = dd.concat([dt[idxb0].Mzams_0,dt[idxb1].Mzams_1])
MergingBHzams = dd.concat([dt[idxm0].Mzams_0,dt[idxm1].Mzams_1])

#Filter and join initial semimajor axis
AllBHa = dd.concat([dt[idx0].a,dt[idx1].a])
BoundBHa = dd.concat([dt[idxb0].a,dt[idxb1].a])
MergingBHa = dd.concat([dt[idxm0].a,dt[idxm1].a])


#Plot
plt.figure(figsize=(10,5))

plt.subplot(1,2,1)
plt.scatter(AllBHzams,AllBH,zorder=1,edgecolor="k",s=30,label="All")
plt.scatter(BoundBHzams,BoundBH,zorder=2,edgecolor="k",s=30, label="Bound")
plt.scatter(MergingBHzams,MergingBH,zorder=3,edgecolor="k",s=30, label="Merging")
plt.plot(np.linspace(0,140),np.linspace(0,140),ls="dashed",c="gray")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("BH mass [M$_\odot$]",fontsize=18)
plt.xlabel("$M\mathrm{zams}$  [M$_\odot$]",fontsize=18)
plt.gca().tick_params(axis='both', which='major', labelsize=18)
plt.legend(fontsize=16)

plt.subplot(1,2,2)
plt.scatter(AllBHa,AllBH,zorder=1,edgecolor="k",s=30,label="All")
plt.scatter(BoundBHa,BoundBH,zorder=2,edgecolor="k",s=30,label="Bound")
plt.scatter(MergingBHa,MergingBH,zorder=3,edgecolor="k",s=30,label="Merging")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Semimajor initial  [R$_\odot$]",fontsize=18)
plt.ylabel("BH mass [M$_\odot$]",fontsize=18)
plt.gca().tick_params(axis='both', which='major', labelsize=18)

plt.tight_layout()
plt.savefig("analysis3.png")
plt.show()
