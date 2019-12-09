This script updates the force field parameters in the tables of the .dms structure file of either receptor or ligand from new force field parameters obtained in text format.

Requirements:
1. file containing the updated parameters in a .txt format
2. The database file <filename.dms> containing individual tables for the force field parameters

Usage:

python update_opls3_param.py <paramfile> <jobname>

jobname: name of the file without the extension

Example:

The example directory has two files:
1. fxr26.txt - contains the new force field parameters
2. fxr26.dms - structure file whose parameters needs to be updated.

The script keeps an untouched copy of the structure file and updates the second copy with new parameters.
The output filename has the naming format <jobname_opls3.dms>

Run the script from the "example" directory:

cd example
python ../update_opls3_param.py fxr26.txt fxr26



