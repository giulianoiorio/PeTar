"""
This example shows how to use RegEX to collect information in the logfile,
in particular it shows how to get the module of the SN kicks for each SN explosion.
The output files are stored in ../sevn_outputs_example/sevn_output_single
"""
import re
import numpy as np

"""1- Create the regex string"""
regex_str="SN;.*?;(?:.*?:){6}(.*?):"
#The string will look for the rows in the logfile regarding the even SN, and will
#get the 7th number in the info part


"""2- Read the logfile and apply regex find all"""
logfile_path="../sevn_outputs_example/sevn_output_single/logfile_0.dat"
with open(logfile_path) as fo:
    ma=re.findall(regex_str,fo.read())
#ma is a list of all the matched output (they are all strings)
print("List of Vkicks taken from re findall")
print(ma)

"""3- Transform a list of string to a numpy array of float"""
vkick=np.array(ma,dtype=float)

"""4- Plot the Vkick distribution"""
import matplotlib.pyplot as plt
plt.hist(vkick,range=(0,600),bins=20)
plt.xlabel("SN Vkick",fontsize=18)
plt.ylabel("$N$",fontsize=18)
plt.yscale("log")
plt.show()
