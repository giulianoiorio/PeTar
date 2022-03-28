Compile SEVN if needed (./compile.sh, Insert the correct SEVN path inside compile.sh)
Insert the correct SEVN path inside run.sh
Run SEVN with ./run.sh 

You will notice that SEVN fails returning an error. 
This is because the input file contains a last addition column with the random seeds. 
In order to enable the rseed input you have to set the parameter RSEED="true" inside the run script. 

Modify the run script  run SEVN again
