"""
This example shows how to use RegEX to collect information in the logfile,
in particular it shows how to get the modulethe number of CE for each events
and add this information to all systems ending their evolution as binary of compact objects
The output files are stored in ../sevn_outputs_example/sevn_output_single
"""
import re
import numpy as np
import pandas as pd  #pandas is already shipped with python anaconda, otherwise install it with pip install pandas


"""1- Count CE for each binary from logfile"""
regex_str=r'B;\d+;(\d+);CE;'

logfile_path="../sevn_outputs_example/sevn_output_single/logfile_0.dat"
with open(logfile_path,"r") as f:
        #Find al the CE occurence and save the ID
        ma = re.findall(regex_str,f.read())
#Create a dataframe with all the ID obtained in the re findall
BID = pd.DataFrame({'ID':np.asarray(ma,dtype=int)})
#Group by ID and count occurences producing a dataframe with the column ID and the column NCE containing the number of CEs
BID= BID.groupby(['ID']).size().reset_index(name='NCE')
print("Dataframe with Number of CE for each Binary (identify by ID) with at least 1 CE")
print(BID)

"""2- Load the output file"""
output_path="../sevn_outputs_example/sevn_output_single/output_0.csv"
df = pd.read_csv(output_path) #df is a pandas dataframe
##Notice if SEVN has been run with the option omode ascii, the fill is called output_0.dat
##and can be load with:
#df = pd.read_csv(output_path,  sep=r"\s+") //read file
#df.rename(columns={'#ID':ID})  //rename ID column to remove the #
print("Loaded file in dataframe")
print(df.head()) #Use head to print just the first few rows


"""3- Take only binaries of compact remnants"""
#Dataframe can be filtered using boolean operations on their columns as done with nunpy array,
#For example, let assume we want to retrieve all the BH-BH bound binaries.
#Remember that the BHs hav RemnantType=6 and bounded binaries have not nan Semimajor
idx_compact_primary= (df.RemnantType_0==4) | (df.RemnantType_0==5) | (df.RemnantType_0==6)
idx_compact_secondary= (df.RemnantType_1==4) | (df.RemnantType_1==5) | (df.RemnantType_1==6)
idxcomp=(idx_compact_primary) & (idx_compact_secondary) & (df.Semimajor.notnull())
#now use the above index to extract a new dataframe of bounde compact binaries
dfcomp=df[idxcomp]
print(f"We found {len(dfcomp)} compact binaries")
print("This is the dataframe ")
print(dfcomp.head())

"""4- Merge the two dataframes"""
dfmerged=pd.merge(dfcomp,BID,on='ID',how="left")
#Notice, we use how="left", because we want that all the binary in the dfcomp dataframe
#are included in the final dataframe, the one with 0 CE events will have a NaN in the NCE column.
#We will set this NaN to 0 in the following rows
dfmerged['NCE'] = dfmerged['NCE'].fillna(0)
dfmerged = dfmerged.astype({"NCE": int})
print("Final merged dataframe")
print(dfmerged.head())

"""5- Plot number of CE"""
import matplotlib.pyplot as plt
bins=(-0.5,0.5,1.5,2.5,3.5)
plt.hist(dfmerged['NCE'],bins=bins)
plt.xlabel("NCE",fontsize=18)
plt.ylabel("$N$",fontsize=18)
plt.yscale("log")
plt.show()
