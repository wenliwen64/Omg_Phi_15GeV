# This script is used to get the phi meson histogram files list only. 
# If the files spread in multiple directories, you have to run this code
# for each directory containing the histogram files.
#
# Example: >>python makefilelist.py 0913_newaa_phi_data/ total_list
# The 0913_newaa_phi_data/ is the directory containing the histogram files.
# The total_list is the files which stores the output result.

import glob 
import os
import string
from sys import argv
from os.path import join

script, input, output = argv

mypath = "/star/u/lwen1990/ucla/v0_19GeV/realoutput/"
file_expression = "*.phi.histo.root"

filepath = mypath + input + file_expression
print filepath
histo_files = glob.glob(filepath);
#histo_files = glob.glob(mypath);
print histo_files

filelist_output = open(output, 'a');

for file in histo_files:
    filelist_output.write(file); 
    filelist_output.write('\n'); 

filelist_output.close();
