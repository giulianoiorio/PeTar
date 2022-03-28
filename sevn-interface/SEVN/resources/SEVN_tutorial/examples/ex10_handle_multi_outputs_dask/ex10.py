"""
This example shows how to use pandas to join the evolve and output files from a SEVN run.
The output files are stored in ../sevn_outputs_example/sevn_output_single
"""
import dask.dataframe as dd
import pandas as pd  #pandas is already shipped with python anaconda, otherwise install it with pip install pandas

"""1- Load the output csv file"""
output_path="../sevn_outputs_example/sevn_output_single/output_*.csv"
dt = dd.read_csv(output_path) #dt is a dask dataframe
##Notice if SEVN has been run with the option omode ascii, the fill is called output_0.dat
##and can be load with:
#dt = dd.read_csv(output_path,  sep=r"\s+") //read file
#dt=df.rename(columns={'#ID':ID})  #rename ID column to remove the #


"""2- Filter the data"""
#Dataframe can be filtered using boolean operations on their columns as done with nunpy array,
#For example, let assume we want to retrieve all the BH-BH bound binaries.
#Remember that the BHs hav RemnantType=6 and bounded binaries have not nan Semimajor
idxBHBH=(dt.RemnantType_0==6) & (dt.RemnantType_1==6) & (dt.Semimajor.notnull())
#now use the above index to extract a new dataframe of bounde BHBH binaries
dtbhbh=dt[idxBHBH]


"""3- Load evolved files """
evolved_path="../sevn_outputs_example/sevn_output_single/evolved_*.dat"
dte = dd.read_csv(evolved_path,  sep=r"\s+") #dt is a dask dataframe
dte = dte.rename(columns={"#ID":"ID"})  #rename ID column to remove the #
#It could be useful to rename the columns name to avoid to have same names of the output files
dte = dte.rename(columns={"Mass_0":"Mzams_0", "Mass_1":"Mzams_1","a":"Semimajor_ini","e":"Eccentricity_ini"})

"""4- Join the filtered dataset with the evolved file"""
merged_dt = dtbhbh.merge(dte,on=["ID","name"],how="inner",suffixes=("","_"))


"""5- Comput all the dask operation and create a pandas dataframe"""
merged_df = merged_dt.compute()

"""6- Save file"""
merged_df.to_csv("results.csv",index=False)

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
