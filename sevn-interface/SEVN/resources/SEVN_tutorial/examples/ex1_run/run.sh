#!/usr/bin/env bash
cmdargument=$1

RED="\033[1;31m"
GREEN="\033[1;32m"
NC="\033[0m" # No Color






#########PARAMETERS#############################

#-------------------------------
#SEVNPATH
#-------------------------------
SEVN="/Users/iogiul/SEVN_dir/SEVN" #Complete path to the SEVN folder
SEVNEXE="${SEVN}/build/exe/sevnB.x" #Complete path to SEVN BSE executable

#-------------------------------
#RUN OPTIONS
#-------------------------------
NTHREADS="1" #Number of OpenMP threads (1 means no parallel threads, sequential execution)
NCHUNK="1" #Evolve Nchunk at time
DTOUT="list"      #If list use the dtout reported in the input list, otherwise use this value for all the stars and binaries (Can be a number in Myr (e.g. 10), a colon separated sequence in Myr (e.g. 10:100:10 goes from 10 Myr to 100 Myr with 10 Myr step, or a list of comma separated numbers in Myr inside curly brackets (e.g. {10,50,40}))


#-------------------------------
#Input and Tables
#-------------------------------
if [ "$cmdargument" = "" ]; then
  LISTBIN="SEVNIC_BSE.in" #Complete path to input file (list of binaries or single stars)
else
  LISTBIN=$cmdargument
fi
IBMODE="new" #Input file format for binaries [new*] [legacy] [sevn1]
TABLES="${SEVN}/tables/SEVNtracks_parsec_AGB" #Complete path to look-up tables
TABLESHE="${SEVN}/tables/SEVNtracks_parsec_pureHe36" #Complete path to look-up tables for pure-He stars
TEND="list"
TSTART="list"
RSEED="false"




#-------------------------------
#OUTPUT
#-------------------------------
OUTPATH="sevn_output/" #Complete path to the output folder (the folder will be automatically created or cleaned if it already exists)
OMODE="csv" #Format for output files [h5] [ascii*] [csv]
NAMEPREX=""  #prefix to add to the name of the systems
LOGLEVEL="critical"  #Log output level: [debug] [info] [warning] [error]
LITPHASES="true" #Use literal phases instead of numbers in output [true] - [false*]
LOGFILE="false" #If true produce the logfile output  [true] - [false*]
SCOL="Worldtime:Mass:Phase:RemnantType" #Additional columns to print in the output file for single stellar evolution runs. Default is empty, but any property of single stars can be added (check names in the Property class)
BCOL="Semimajor:Eccentricity:BEvent" #Additional columns to print in the output file for binary stellar evolution runs. Default is empty, but any property of binary stars can be added (check names in the BinaryProperty class)


#-------------------------------
#METALLICITY
#-------------------------------
Z="list"  #Stellar metallicity - [*list][number]. If list use the Z in the input file otherwise overwrite all the Zs.

#-------------------------------
#PRESCRIPTIONS
#-------------------------------
WINDSMODE="hurley" #Prescriptions for wind accretion and the associated orbital changes - [hurley*] [disabled]
TIDES="tides_simple" #Prescriptions for tides - [simple*] [disabled]
GWMODE="peters" #Prescriptions for gravitational-wave decay - [peters*] [disabled]
RLMODE="hurley_rl" #Prescriptions for Roche-Lobe overflow and mass transfer/accretion - [hurley_rl*] [disabled]
CEMODE="energy" #Prescriptions for common-envelope evolution - [hurley*] [disabled]
MIXING="simple" #Prescriptions to for mixing - [simple*] [disabled]
SNORBCHANGE="hurley" #Prescriptions for orbital changes after SN kicks - [hurley*] [disabled]
SNKICKS="unified" #Prescriptions for SN kicks - [unified*] [hobbs] [zeros]
SUPERNOVA="list" #Prescription for the SN explosion mechanism - [list*] [rapid] [delayed] [compactness]
#---------OPTIONS---------------
SQHE="false" #If true enable the Quasi Homogeneous Evolution  after a RLO mass transfer following Elrdige&Stanway11
TABCONV="true" #If true estimate the properties of the convective envelope using the tables (xxxconv.dat)
TABXSUP="false"  #If true estimate the superficial abundance  using the tables (xxxsup.dat)
TABINERTIA="true" #If true estimate the properties of the stellar inerita using the tables (inertia.dat)
TABRHE="true" #If true estimate the properties of the HE core radius using the tables (rhe.dat)
TABRCO="true" #If true estimate the properties of the CO core radius using the tables (rhe.dat)

#-------------------------------
#Parameters
#-------------------------------
#-------GW-------#
GWTSHOLD="1" #Enable GW decay if GW_time_decay < GWTSHOLD*Hubble_time
#-------RLO-------#
RLOEDD="1" #Eddington factor to limit accretion on a compact object (>1 means super-Eddington)
RLOEPSNOVA="0.001" #Fraction of accreted matter retained in nova eruption
RLOMACCR="0.7" #Fraction of the mass lost by the primary that is accreted onto the secondary during RLO
RLOGAM="-1" #Angular momentum lost during RLO. [-1]: from the primary, [-2]: from the secondary, [>0]: from the system
RLOCIRC="false" #If true the orbit is circularised at the onset of the RLO
RLOQCRIT="Hurley_Webbink"
#-------SN-------#
MCHANDRA="1.41" #Chandrasekar mass limit for WD formation
SNLOWECSN="1.38" #Minimum value for the CO mass to go ECSN
SNLOW="1.44" #Minimum CO value for the CO mass to go SN (i.e. max CO mass for ECSN)
SNC25TS="-1" #csi25 parameter threshold for explosion/implosion decision
SNCOMPFB="0.9" #Fallback fraction for explosion in the compact SN option
#-------WINDS-------#
WALPHA="1.5" #alpha factor to tune the amount of wind accretion (Eq.6 Hurley+02)
WBETA="0.125" #beta factor to tune wind velocity (Eq.9 in Hurley+02)
#-------CE-------#
CEALPHA="5" #alpha in binding energy (Eq. 73 in Hurley02)
CELAM="-1" #if >0 Constant Lambda in binding energy (Eq. 69 in Hurley02). If -1 use Lambda from Claeys et al. 2014 (Appendix A).
CELAMFTH="1"  #Fraction of internal energy that goes to the binding energy. (equivalent to lambda_ion in Claeys et al. 2014 (Appendix A)). Used only if star_lambda=-1.
CEKCE="1"  #Fraction of non core mass  participating to the CE (e.g. envelope of giants) retained after the CE coalescence.If -1, use a rescaled version of eq. 77 In Hurley
CEKNCE="1" #Fraction of non core mass not participating to the CE (e.g. a MS star) retained after the CE coalescence.If -1, use the eq. 77 in Hurley 2002 (ce_kce is ignored)
#-------NS-------#
NSMAX="3.0" #Maximum NS mass
NSMAGTSCALE="1000" #Magnetic field decay timescale in Myr
NSMAGMSCALE="0.15" #Magnetic field decay mass-scale in Msun

#ev_naked_tshold:

# -------------------------------
#Parameters for jumping onto new tracks
#-------------------------------
JTEMAX="0.005" #Maximum relative error in mass when jumping on a new track
JTMAXDM="1.2" #Maximum new ZAMS tested when jumping on a new track (Mzams_new_max = Mzams_old + JTMAXDM*DM_accreted_or_donated)
JTMINDM="0" #Minimum new ZAMS tested when jumping on a new track (Mzams_new_min = Mzams_old + JTMINDM*DM_accreted_or_donated)
JTMAXITER="10" #Maximum number of iterations for reaching convergence
JTDMTSHOLD="0.01" #Maximum relative change in total mass for not changing track


# -------------------------------
#Other parameters
#-------------------------------
MAXREP="50" #Maximum number of repetitions allowed in the sse and bse. If we reach this number an error is raised
NAKEDTS="1E-4" #Mass difference threshold (Msun) between envelope and core to set a star as nakedHe or nakedCO.
INITERRSTOP="false" #If true terminate the run when a error on a system initialisation is thrown
SMAXCO="false" #If true the first time a star develops a CO core, we set the maximum CO core Mass for SSE as the last value of the interpolating tracks
SMINHE="false" #If true the first time a star develops a CO core, we set the minimum HE core Mass for SSE as the last value of the interpolating tracks

PSTRING="-nthreads $NTHREADS \
-list $LISTBIN \
-ibmode $IBMODE \
-tables $TABLES \
-tables_HE $TABLESHE \
-omode $OMODE \
-o $OUTPATH \
-wmode $WINDSMODE \
-tmode $TIDES \
-gwmode $GWMODE \
-rlmode $RLMODE \
-cemode $CEMODE \
-mixmode $MIXING \
-kmode $SNORBCHANGE \
-sn_kicks $SNKICKS \
-snmode $SUPERNOVA \
-gw_tshold $GWTSHOLD \
-jtrack_h_err_rel_max $JTEMAX \
-jtrack_max_dm_factor $JTMAXDM \
-jtrack_max_iteration $JTMAXITER \
-jtrack_min_dm_factor $JTMINDM \
-jtrack_tshold_dm_rel $JTDMTSHOLD \
-rlo_eddington_factor $RLOEDD \
-rlo_eps_nova $RLOEPSNOVA \
-rlo_f_mass_accreted $RLOMACCR \
-rlo_gamma_angmom $RLOGAM \
-rlo_qcrit $RLOQCRIT \
-sn_Mchandra $MCHANDRA \
-sn_co_lower_ecsn $SNLOWECSN \
-sn_co_lower_sn $SNLOW \
-sn_max_ns_mass $NSMAX \
-w_alpha $WALPHA \
-w_beta $WBETA \
-io_literal_phases $LITPHASES \
-log_level $LOGLEVEL \
-dtout $DTOUT \
-ev_Nchunk $NCHUNK \
-ev_max_repetitions $MAXREP \
-tf $TEND \
-tini $TSTART \
-rseed $RSEED \
-rlo_circularise $RLOCIRC \
-rlo_QHE $SQHE \
-tabuse_envconv $TABCONV \
-tabuse_Xsup $TABXSUP \
-tabuse_inertia $TABINERTIA \
-tabuse_rco  $TABRCO \
-tabuse_rhe  $TABRCO \
-ce_alpha $CEALPHA \
-star_lambda $CELAM \
-star_lambda_fth $CELAMFTH \
-ce_kce $CEKCE \
-ce_knce $CEKNCE \
-ns_magnetic_tscale $NSMAGTSCALE \
-ns_magnetic_mscale $NSMAGMSCALE \
-sn_compact_csi25_tshold $SNC25TS \
-sn_compact_fallback $SNCOMPFB \
-ev_naked_tshold $NAKEDTS \
-Z $Z  \
-io_logfile $LOGFILE \
-initerror_stop $INITERRSTOP \
-ev_set_maxCO $SMAXCO \
-ev_set_minHE $SMINHE"


#Check SCOL, BCOL and NAMEPREX, if they are empty do not add it, 
#otherwise they will generate errors (SEVN expect a non-empty value after an option)
#the [[:blank:]] is a command to remove all the white spaces.

if [ "${SCOL//[[:blank:]]/}" != ""   ]; then
  PSTRING="$PSTRING -scol $SCOL"
fi

if [ "${BCOL//[[:blank:]]/}" != "" ]; then
  PSTRING="$PSTRING -bcol $BCOL"
fi

if [ "${NAMEPREX//[[:blank:]]/}" != "" ]; then
   PSTRING="$PSTRING -name_prefix $NAMEPREX"
fi



LAUNCH="$SEVNEXE $PSTRING"
DATETIME=$(date '+%^B %d, %Y @ %H:%M:%S');

mkdir -p $OUTPATH


#This is where we have the actual launch command
echo -e "${GREEN}SEVN is running...${NC}"
eval $LAUNCH
retVal=$?

if [ $retVal -eq 0 ]; then
    echo -e "${GREEN}SEVN execution was SUCCESSFUL!"
    echo -e "The output can be found in $OUTPATH${NC}"
else
    echo -e "${RED}SEVN execution has FAILED!${NC}"
    exit 1
fi

echo "$DATETIME $HOSTNAME [SEVN launch string]: $LAUNCH" > $OUTPATH/launch_line.txt
