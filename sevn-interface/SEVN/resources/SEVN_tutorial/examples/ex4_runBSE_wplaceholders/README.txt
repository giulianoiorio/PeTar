Compile SEVN if needed (./compile.sh, Insert the correct SEVN path inside compile.sh)
Insert the correct SEVN path inside run.sh

You will notice that SEVN fails returning an error. 
This is because in the initial conditions (SEVNIC_BSE_placeholder.in) there are placeholders for 
Z, SN, start, tend, dtout.

Modify the run script to replace the placeholders and run SEVN again
