# Module to setup the Sofia-2 pipeline
# Generates Sofia parameter files with input subregions obtained by slicing the
# FITS cube along freq channels
# Based on s2p_setup  https://github.com/SoFiA-Admin/s2p_setup

import sys
import os
import math
import configparser
from astropy.io import fits


# Function to search for substring in list
def substr_search(input_list, substr):
    for i, s in enumerate(input_list):
        if substr in s:
            return i
    return -1

def gen_par_files(template_file, output_db_name):

# Read settings from configuration file
    config = configparser.ConfigParser();
    success = config.read("pipeline_setup.ini");

    if(len(success) == 0):
        sys.stderr.write("Error: Failed to read config file: pipeline_setup.ini\n");
        sys.exit(1);

#TODO Assumption that entire data cube needs to be processed, Add boundary configuration

    number_of_slices = int(config["region"]["number_of_slices"])
    overlap_spec     = int(config["region"]["overlap_spec"])


#TODO Read from FITS, Dummy values for testing

    x_min = 0
    y_min = 0
    z_min = 0
    x_max = 500
    y_max = 500
    z_max = 500

#TODO Add RAM and region size calculations

    slice_size_z = int(math.floor(float(z_max - z_min + 1 +(number_of_slices - 1)*overlap_spec) / float(number_of_slices)))

# Read template parameter file
    try:
        with open(template_file) as par_file:
            template_par = par_file.readlines()
    except:
        sys.stderr.write("Error: Failed to read template parameter file sofia.par.\n")
        sys.exit(16)


# Create set of parameter files

    for z in range(number_of_slices):

        par = template_par[:]

        z1 = z_min + int(math.floor(z * (slice_size_z - overlap_spec)))
        z2 = z1 + slice_size_z

        if(z1 < 0): z1 = 0
        if(z2 > z_max or z == number_of_slices): z2 = z_max

        i = substr_search(par, "input.region")

        if(i < 0): par.append("input.region  =  {0:d},{1:d},{2:d},{3:d},{4:d},{5:d}\n".format(x_min, x_max, y_min, y_max, z1, z2));
        else: par[i] = "input.region  =  {0:d},{1:d},{2:d},{3:d},{4:d},{5:d}\n".format(x_min, x_max, y_min, y_max, z1, z2);

        i = substr_search(par, "output.filename");
        if(i < 0): par.append("output.filename  =  {0}_{1:03d}\n".format(output_db_name, z));
        else: par[i] = "output.filename  =  {0}_{1:03d}\n".format(output_db_name, z);

        # Dump parameters into new file
        filename = "sofia_{0:03d}.par".format(z);
        sys.stdout.write("  {0}\n".format(filename));
        try:
            with open(filename, "w") as par_file:
                for item in par:
                    par_file.write("{0}".format(item));
        except:
            sys.stderr.write("Error: Failed to write output parameter file: {}\n".format(filename));
            sys.exit(1);
