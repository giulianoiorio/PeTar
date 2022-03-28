"""
This example show how to read and apply  a simple filter to the output.csv from SEVN using pandas.
The output files are store in ../sevn_outputs_example/sevn_output_single
"""

import pandas as pd  #pandas is already shipped with python anaconda, otherwise install it with pip install pandas

"""1- Load the output csv file"""
output_path="../sevn_outputs_example/sevn_output_single/output_0.csv"
df = pd.read_csv(output_path) #df is a pandas dataframe
##Notice if SEVN has been run with the option omode ascii, the fill is called output_0.dat
##and can be load with:
#df = pd.read_csv(output_path,  sep=r"\s+") //read file
#df.rename(columns={'#ID':ID})  //rename ID column to remove the #
print("Loaded file in dataframe")
print(df.head()) #Use head to print just the first few rows
print("Columns in the daframe")
print(df.columns)


"""2- Access data in the dataframe"""
#Data in the daframe can be accessed through the name of column with the dot formalism,
#e.g. df.name_column or with the dictonary formalism df["name_column"]
#For example if we want to  get  Mass_0 vs Mass_1
Mass = df.Mass_0 #or df["Mass_0"]
Semi = df["Mass_1"] #or df.Mass_1
#plot it
import matplotlib.pyplot as plt
plt.scatter(Mass,Semi)
plt.xlabel("Mass_0 [Msun]")
plt.ylabel("Mass_1 [Msun]")
plt.show()

"""3- Filter the data"""
#Dataframe can be filtered using boolean operations on their columns as done with nunpy array,
#For example, let assume we want to retrieve all the BH-BH bound binaries.
#Remember that the BHs hav RemnantType=6 and bounded binaries have not nan Semimajor
idxBHBH=(df.RemnantType_0==6) & (df.RemnantType_1==6) & (df.Semimajor.notnull())
#now use the above index to extract a new dataframe of bounde BHBH binaries
dfbhbh=df[idxBHBH]
print(f"We found {len(dfbhbh)} BHBH binaries")
print("This is the dataframe ")
print(dfbhbh.head())

Mass = dfbhbh.Mass_0 #or df["Mass_0"]
Semi = dfbhbh["Mass_1"] #or df.Mass_1
#plot it
import matplotlib.pyplot as plt
plt.scatter(Mass,Semi)
plt.xlabel("BH_Mass_0 [Msun]")
plt.ylabel("BH_Mass_1 [Msun]")
plt.show()
