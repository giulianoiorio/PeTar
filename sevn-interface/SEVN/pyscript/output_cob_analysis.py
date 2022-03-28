#!/usr/bin/env python
"""
SCRIPT: output_cob_analysis.py
PURPOSE: Take all the output of a sevn simulation filter compact object binaries
and save the following files:
	- BHBH.csv: containing all the bound BHBH binaries
	- BHBHm.csv: containing all the bound BHBH binaries that will merge within an Hubble time (14000 My)
	- NSNS.csv: containing all the bound NSNS binaries
	- NSNSm.csv: containing all the bound NSNS binaries that will merge within an Hubble time (14000 My)
	- BHBH.csv: containing all the bound BHNS binaries
	- BHBHm.csv: containing all the bound BHNSm binaries that will merge within an Hubble time (14000 My)

There are two ways to use this script,
****1- Analyse a single sevn output folder*********
	USAGE;  ./output_cob_analysis.py [-h] [-n NPROC] [-o OUT] input
	e.g. ./output_cob_analysis.py sevn_output -n 2
	positional arguments:
	  input                 path to the sevn folder with the input to analyse

	optional arguments:
	  -h, --help            show this help message and exit
	  -n NPROC, --nproc NPROC
							Number of parallel processes to use [2]
	  -o OUT, --output OUT  Output folder where to store the files, if None use the same input folder [None]
*******************************************************************

****2- Analyse multiple sevn output folders*********
	usage: ./output_cob_analysis.py -m [-h] [-n NPROC]   [--froot FROOT] [--subfolder SUBFOLDER] input
	e.g. ./output_cob_analysis.py -m   -n 2  --froot sevn_output --subfolder 0 main_folder
		In this case, the script will analyse all the data folder starting with  name starting with
		sevn_output inside the main_folder, the output data are inside the subfolder with name 0, e.g.
		sevn_output_Z02/0, the output are automatically saved inside the root sevn output folders
		(e.g. sevn_output_Z02)

	positional arguments:
	  input                 root folder containing all the sevn output folders to analyse

	optional arguments:
	  -h, --help            show this help message and exit
	  -n NPROC, --nproc NPROC
							Number of parallel processes to use [2]
	  -m, --multifolders    If true analyse a list of output folders inside the path specified
	  --froot FROOT         common root for multiple sevn output analysis, used only if option -m is enabled [sevn_output]
	  --subfolder SUBFOLDER
							Subfolder containg the SEVN output, used only if option -m is enabled [None]on input and with folder root specified with
							the froot parameter [False]

REQUIRED PACKAGE:
    -numpy
    -pandas

V 1.0: 04/01/22 Giuliano Iorio: First version
V 1.1: 05/01/22 Giuliano Iorio: Added the option to analyse multiple sevn output folders
"""

import numpy as np
import pandas as pd
import re
import glob
from itertools import repeat
from multiprocessing import Pool
import argparse
import os

class SEVN_output:

	def __init__(self,outfolder):

		self.outfolder = outfolder

	@staticmethod
	def get_history(logfile):
		"""
		Extract information about the Binary event for each system
		:param logfile: Path to the logfile to analyse
		:return: a pandas dataframe with columns: name, ID, Events,
		the Events column contains string with all the keyword of the events separated by a :,
			RB= RLO BEGIN
			RE= RLO END
			C= Common Envelope
			K= Collision
			M= Merger
			W= Swallowed
			S= Supernova Explosion
		"""
		dic_map={"RLO_BEGIN":"RB","RLO_END":"RE","CE":"C","COLLISION":"K",
				 "MERGER":"M","SWALLOWED":"W","BSN":"S"}

		list_map = lambda x: ":".join([dic_map[_x] for _x in x])

		with open(logfile,"r") as fo:
			i=fo.read()
			ma = re.findall("B;([0-9]_?[0-9]*);([0-9]*);(.*);.*;",i)
		na = np.array(ma)
		df = pd.DataFrame(na,columns=["name","ID","Events"])
		df = df.groupby(["ID","name"]).agg(list_map).reset_index()
		return df.astype({"ID":"int64"})

	@staticmethod
	def get_COBs(outfile):
		"""
		Filter all the bound Compact object binaries (BHBH BHNS NSNS) from an outputfile.
		Notice this function assume that the simulation is stopped when both stars are remnant, therefore it gives more than
		one result for systems if the evolution is not stopped at the remnant formation.
		:param outfile: Path to the output file to analyse
		:return: a pandas dataframe containing the all COB binaries
		"""
		df = pd.read_csv(outfile)
		idx_BBH = (df.RemnantType_0 == 6) & (df.RemnantType_1 == 6)
		idx_BNS = ((df.RemnantType_0 == 4) | (df.RemnantType_0 == 5)) & ((df.RemnantType_1 == 4) | (df.RemnantType_1 == 5))
		idx_BNB = ((df.RemnantType_0 == 6) & ((df.RemnantType_1 == 4) | (df.RemnantType_1 == 5))) | (
				(df.RemnantType_1 == 6) & ((df.RemnantType_0 == 4) | (df.RemnantType_0 == 5)))
		idx = (idx_BBH | idx_BNS | idx_BNB) & (df.Semimajor > 0)
		dff = df[idx]

		return dff.astype({"name":object})

	@staticmethod
	def get_IC(outfile):
		"""

		:param outfile: Path to the evolved file to analyse
		:return: a pandas dataframe with the column:
			- ID: system id
			- name: sytem name
			- Mzams_0: initial zams mass of the first star
			- Mzams_1: initial zams mass of the second star
			- a: initial semimajor axis
			- e: initial eccentricity
			- Z: metallicity
		"""
		df = pd.read_csv(outfile, sep="\s+")
		df = df.rename(columns={"#ID": "ID",
								"Mass_0": "Mzams_0",
								"Mass_1": "Mzams_1",
								"a": "Semimajor_ini",
								"e": "Eccentricity_ini",
								"Z_0": "Z"})

		return df

	@staticmethod
	def make_COB_analysis_single(file_id, outpath):
		"""
		Return a pandas dataframe containg all the bound compact object binaries with information from the output, evolved and logfile
		:param file_id: number of thread that generate the output to analyse
		:param outpath: Path to the folder containing the sevn output
		:return:
		"""
		logfile = f"{outpath}/logfile_{file_id}.dat"
		outfile = f"{outpath}/output_{file_id}.csv"
		evolvedfile = f"{outpath}/evolved_{file_id}.dat"
		# History
		dfh = SEVN_output.get_history(logfile)
		# Evolution
		dfe = SEVN_output.get_COBs(outfile)

		# Filter columns
		cols_in_file=dfe.columns
		wanted_columns=["ID", "name", "BWorldtime", "Mass_0", "Radius_0", "Zams_0", "Phase_0", "RemnantType_0", "Mass_1", "Radius_1",
						"Phase_1", "Zams_1", "RemnantType_1", "Semimajor", "Eccentricity", "GWtime"]
		effective_columns = [col for col in wanted_columns if col in cols_in_file]
		dfe=dfe[effective_columns]
		dfe["name"] = dfe["name"].astype('str')
		# IC
		dfi = SEVN_output.get_IC(evolvedfile)[["ID", "name", "Mzams_0", "Mzams_1", "Semimajor_ini", "Eccentricity_ini", "Z"]]
		dfi["name"] = dfi["name"].astype('str')

		dff = dfe.merge(dfh, on=["ID", "name"], how="left")
		dff = dff.merge(dfi, on=["ID", "name"], how="left")

		return dff

	def make_COB_analysis(self,nproc=4):
		"""

		:param nproc: number of parallel process to use
		:return:
		"""
		files = glob.glob(self.outfolder + "/logfile_*")
		idlist = [int(x.split("logfile_")[1].split(".")[0]) for x in files]

		with Pool(nproc) as p:
			dfl = p.starmap(SEVN_output.make_COB_analysis_single, zip(idlist, repeat(self.outfolder)))

		return pd.concat(dfl)

	@staticmethod
	def split_COBs(df):

		idxbbh = (df.RemnantType_0==6) & (df.RemnantType_1==6) & (df.Semimajor.notnull())
		idxbns = ( (df.RemnantType_0==4) | (df.RemnantType_0==5) ) &  ((df.RemnantType_1==4) | (df.RemnantType_1==5) ) & (df.Semimajor.notnull())
		idxbhns = ( ( ((df.RemnantType_0==4) | (df.RemnantType_0==5)) & (df.RemnantType_1==6)  ) \
					| ( ((df.RemnantType_1==4) | (df.RemnantType_1==5)) & (df.RemnantType_0==6)  )  ) & (df.Semimajor.notnull())

		dfbbh  = df[idxbbh]
		dfbns  = df[idxbns]
		dfbhns = df[idxbhns]

		return dfbbh, dfbns, dfbhns

	@staticmethod
	def merging(df,tshold=14000):

		idxmerging = df.BWorldtime + df.GWtime < tshold

		return df[idxmerging]



def analyse_sevn_output(input_folder, nproc=2, output_older=None):
	"""
	Analyse the file of a sevn output filtering the compact objectes binaries mixing information from output, evolved an logfiles.
	6 csv files will be saved:
	- BHBH.csv: containing all the bound BHBH binaries
	- BHBHm.csv: containing all the bound BHBH binaries that will merge within an Hubble time (14000 My)
	- NSNS.csv: containing all the bound NSNS binaries
	- NSNSm.csv: containing all the bound NSNS binaries that will merge within an Hubble time (14000 My)
	- BHBH.csv: containing all the bound BHNS binaries
	- BHBHm.csv: containing all the bound BHNSm binaries that will merge within an Hubble time (14000 My)
	:param input_folder: Name of sevn output folder containt the output, evolved and logfiles
	:param nproc: number of parallel processes to use
	:param output_older: name of the folder where to store the output file, if None the same input folder will be used
	:return: 0
	"""
	if output_older is None:
		output_older = input_folder

	if not os.path.isdir(output_older):
		os.makedirs(output_older)

	so = SEVN_output(input_folder)
	df = so.make_COB_analysis(nproc)
	dfbbh, dfbns, dfbhns = so.split_COBs(df)

	dfbbh.to_csv(output_older + "/BHBH.csv", index=False)
	dfbhbhm = so.merging(dfbbh)
	dfbhbhm.to_csv(output_older+ "/BHBHm.csv", index=False)

	dfbns.to_csv(output_older + "/NSNS.csv", index=False)
	dfbnsm = so.merging(dfbns)
	dfbnsm.to_csv(output_older + "/NSNSm.csv", index=False)

	dfbhns.to_csv(output_older + "/BHNS.csv", index=False)
	dfbhnsm = so.merging(dfbhns)
	dfbhnsm.to_csv(output_older + "/BHNSm.csv", index=False)

	return 0

def analyse_multi_sevn_output(main_folder,  nproc=2, sevn_root="sevn_output", subfolder_name=None):
	"""
	Analyse the file of multiple sevn outputs folder having the same root name.
	For each simulation folder the output files are analysed filtering the compact objectes binaries mixing information from output, evolved an logfiles.
	6 csv files will be saved:
	- BHBH.csv: containing all the bound BHBH binaries
	- BHBHm.csv: containing all the bound BHBH binaries that will merge within an Hubble time (14000 My)
	- NSNS.csv: containing all the bound NSNS binaries
	- NSNSm.csv: containing all the bound NSNS binaries that will merge within an Hubble time (14000 My)
	- BHBH.csv: containing all the bound BHNS binaries
	- BHBHm.csv: containing all the bound BHNSm binaries that will merge within an Hubble time (14000 My)
	:param main_folder: path to the folder containing the multiple sevn output folders
	:param nproc: number of parallel processes to use
	:param sevn_root: common root name of the sevn output folders, e.g. if we have sevn_out1 sevn_out2,
	the root name should be sevn_out
	:param subfolder_name: Sometime the outputs are not contained directly in the main sevn output folder but thera a number
	of subfolders to consider, this is the name of the subfolder we want to analyse.
	:return: 0
	"""

	folders = glob.glob(main_folder+"/"+sevn_root+"*")

	for folder in folders:
		print("Analysing folder ",folder, flush=True)
		fname=folder
		if subfolder_name is not None:  fname = folder+"/"+subfolder_name
		try:
			analyse_sevn_output(fname, nproc=nproc, output_older=folder)
			print("Done with folder ",folder, flush=True)
		except Exception as e:
			print("Error analyse folder ",folder," with message: ",file=sys.stderr)
			print(e,file=sys.stderr)


	return 0

if __name__=="__main__":

	info = "	Analyse the file of a sevn output filtering the compact objectes binaries mixing information from output, evolved an logfiles.\n"
	info += "6 csv files will be saved:\n"
	info += "- BHBH.csv: containing all the bound BHBH binaries\n"
	info += "- BHBHm.csv: containing all the bound BHBH binaries that will merge within an Hubble time (14000 My)\n"
	info += "- NSNS.csv: containing all the bound NSNS binaries\n"
	info += "- NSNSm.csv: containing all the bound NSNS binaries that will merge within an Hubble time (14000 My)\n"
	info += "- BHNS.csv: containing all the bound BHNS binaries\n"
	info += "- BHNSm.csv: containing all the bound BHNS binaries that will merge within an Hubble time (14000 My)"


	parser = argparse.ArgumentParser(description=info)
	parser.add_argument('input', type=str, help='path to the sevn folder with the input to analyse, if option -m is enabled, this'
												'is the root folder containing all the sevn output folders to analyse')
	parser.add_argument('-n','--nproc', dest='nproc', default=2, type=int, help="Number of parallel processes to use [2]")
	parser.add_argument('-o', '--output', dest='out', default=None, type=str, help="Output folder where to store the files, if None the input folder is used [None]")
	parser.add_argument('-m', '--multifolders', dest='multifolders', action='store_true', default=False, help="If true analyse a list of output folders inside the path specified")
	parser.add_argument('--froot', dest='froot', default="sevn_output", type=str, help="common root for multiple sevn output analysis, used only if option -m is enabled [sevn_output]")
	parser.add_argument('--subfolder', dest='subfolder', default=None, type=str, help="Subfolder containg the SEVN output, used only if option -m is enabled [None]"
																					  "on input and with folder root specified with the froot parameter [False]")
	args = parser.parse_args()

	if (args.multifolders):
		analyse_multi_sevn_output(args.input,args.nproc,args.froot,args.subfolder)
	else:
		analyse_sevn_output(args.input,args.nproc,args.out)


