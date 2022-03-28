"""
This example shows how to use pandas to join the evolve and output files from a SEVN run.
The output files are stored in ../sevn_outputs_example/sevn_output_single
"""

import pandas as pd  #pandas is already shipped with python anaconda, otherwise install it with pip install pandas

"""1- Load the output csv file"""
output_path="../sevn_outputs_example/sevn_output_single/output_0.csv"
df = pd.read_csv(output_path) #df is a pandas dataframe
##Notice if SEVN has been run with the option omode ascii, the fill is called output_0.dat
##and can be load with:
#df = pd.read_csv(output_path,  sep=r"\s+") //read file
#df=df.rename(columns={'#ID':ID})  #rename ID column to remove the #
print("Loaded file in dataframe")
print(df.head()) #Use head to print just the first few rows
print("Columns in the daframe")
print(df.columns)

"""2- Filter the data"""
#Dataframe can be filtered using boolean operations on their columns as done with nunpy array,
#For example, let assume we want to retrieve all the BH-BH bound binaries.
#Remember that the BHs hav RemnantType=6 and bounded binaries have not nan Semimajor
idxBHBH=(df.RemnantType_0==6) & (df.RemnantType_1==6) & (df.Semimajor.notnull())
#now use the above index to extract a new dataframe of bounde BHBH binaries
dfbhbh=df[idxBHBH]
print(f"We found {len(dfbhbh)} BHBH binaries")
print("This is the dataframe ")
print(dfbhbh.head())

"""3- Load evolved files """
evolved_path="../sevn_outputs_example/sevn_output_single/evolved_0.dat"
dfe = pd.read_csv(evolved_path,  sep=r"\s+") #df is a pandas dataframe
dfe=dfe.rename(columns={"#ID":"ID"})  #rename ID column to remove the #
#It could be useful to rename the columns name to avoid to have same names of the output files
dfe = dfe.rename(columns={"Mass_0":"Mzams_0", "Mass_1":"Mzams_1","a":"Semimajor_ini","e":"Eccentricity_ini"})
print("Loaded evolved file in dataframe")
print(dfe.head()) #Use head to print just the first few rows
print("Columns in the daframe")
print(dfe.columns)




"""4- Join the filtered dataset with the evolved file"""
merged_df = dfbhbh.merge(dfe,on=["ID","name"],how="inner",suffixes=("","_"))
print("Created mergerd eataframe")
print(merged_df.head()) #Use head to print just the first few rows
print("Columns in the daframe")
print(merged_df.columns)

##Now plot the Initial semimajor axis of the system (from the evolved file) vs the
#Semimajor axis at the time of BHBH formation (from the output file)
import matplotlib.pyplot as plt
Semimajor_ini=merged_df.Semimajor_ini
Semimajor_end=merged_df.Semimajor
plt.scatter(Semimajor_ini,Semimajor_end)
plt.xlabel("Semimajor axis (t=0) [Rsun]")
plt.ylabel("Semimajor axis at BHBH formation [Rsun]")
plt.xscale("log")
plt.yscale("log")
plt.show()
