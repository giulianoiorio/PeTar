import dask.dataframe as dd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#Load file
dt=dd.read_csv("sevn_output/output_*.csv")
#Give a look to the columns
print(dt.columns)

time_ranges=((0,1),(1,10),(10,50),(300,500))

plt.figure(figsize=(12,10))

for i,timer in enumerate(time_ranges):
    tlow,tup=timer

    plt.subplot(2,2,i+1)
    #Filter only stars
    idxt=(dt.BWorldtime>=tlow) & (dt.BWorldtime<=tup)
    idx0 = (dt.Phase_0<7) & idxt
    idx1 = (dt.Phase_1<7) & idxt
    Luminosity = dd.concat([dt.Luminosity_0[idx0],dt.Luminosity_1[idx1]])
    Temperature = dd.concat([dt.Temperature_0[idx0],dt.Temperature_1[idx1]])
    plt.hexbin(np.log10(Temperature),np.log10(Luminosity),bins=100,cmap="plasma",mincnt=1)
    cbar=plt.colorbar(pad=0)
    cbar.ax.tick_params(axis='both', which='major', labelsize=16)
    cbar.set_label(label="$N$",size=15)
    plt.xlim(5.5,3.5)
    plt.xlabel("$\log T/\mathrm{K}$",fontsize=18)
    plt.ylabel("$\log L/\mathrm{L_\odot}$",fontsize=18)
    plt.gca().tick_params(axis='both', which='major', labelsize=18)
    plt.ylim(0.5,7)
    plt.title(f"{tlow}<Age/Myr<{tup}",fontsize=20)


plt.tight_layout()
plt.savefig("analysis0.png")
plt.show()
