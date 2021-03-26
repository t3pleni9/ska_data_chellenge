# Helper class to setup the Sofia-2 pipeline
# Based on s2p_setup  https://github.com/SoFiA-Admin/s2p_setup

import sys
import os
import math
import configparser
from astropy.io import fits
from .sofia_params import SofiaParams

class PipelineSetup:
    def __init__(self):
        pass

    # Function to search for substring in list
    def substr_search(self, input_list, substr):
        for i, s in enumerate(input_list):
            if substr in s:
                return i
        return -1

    # Function to generate Sofia parameter files for sub-regions split along x and y axes
    def gen_par_files(self, config):

        #TODO move the config reading check to pipeline init
        # Read settings from configuration file
        # config = configparser.ConfigParser();
        # success = config.read("pipeline_setup.ini");
    
        # if(len(success) == 0):
        #     sys.stderr.write("Error: Failed to read config file: pipeline_setup.ini\n");
        #     sys.exit(1);
    
        #TODO Assumption that entire data cube needs to be processed, Add boundary configuration(if necessary)

        input_file         = config["setup"]["input_file"]
        template_file      = config["setup"]["template_file"]
        par_file_directory = config["setup"]["par_file_directory"]
        output_directory   = config["setup"]["output_directory"]
        output_db_name     = config["setup"]["output_db_name"]

        region_size   = int(config["region"]["region_size"])
        overlap_spat  = int(config["region"]["overlap_spat"])

        # Open FITS file
        try:
            hdu = fits.open(input_file);
        except:
            sys.stderr.write("Error: Failed to open FITS file: {}\n".format(input_file));
            sys.exit(1);
    
        # Extract header information
        header = hdu[0].header;
    
        bitpix = int(header["BITPIX"]);
        word_size = int(abs(bitpix) / 8);
        if(bitpix > 0):
            wordsize = 4; # Assume 32-bit for integer arrays
            sys.stderr.write("Warning: Data cube is of integer type; assuming 32 bit per pixel.\n");
    
        naxis = int(header["NAXIS"]);
        if(naxis < 3 or naxis > 4):
            sys.stderr.write("Error: Data cube is not three-dimensional.\n");
            hdu.close();
            sys.exit(1);
    
        nx = int(header["NAXIS1"]);
        ny = int(header["NAXIS2"]);
        nz = int(header["NAXIS3"]);
    
        if(nz == 1 and naxis == 4):
            nz = int(header["NAXIS4"]);
            sys.stderr.write("Warning: Swapping 3rd and 4th axis of 4D cube.\n");
    
        # Close FITS file
        hdu.close();
    
        # Entire FITS cube
        x_min = 0;
        y_min = 0;
        z_min = 0;
        x_max = nx - 1;
        y_max = ny - 1;
        z_max = nz - 1;
    
        #TODO Add RAM and region size calculations

        # Number of pixels per region
        pixels_per_region = int(math.floor((float(1024 * 1024 * 1024 * region_size) / float(word_size))))

        # Nominal size of each subregion in x and y axes
        size_xy = int(math.floor(float(pixels_per_region)/float(z_max - z_min + 1)) ** (1.0 / 2.0))

        # Number of regions in x and y
        n_reg_x = int(math.ceil(float(x_max - x_min + 1) / float(size_xy)));
        n_reg_y = int(math.ceil(float(y_max - y_min + 1) / float(size_xy)));

        # Size of each region in x and y with overlap
        size_x = int(math.floor(float(x_max - x_min + 1 + (n_reg_x - 1) * overlap_spat) / float(n_reg_x)));
        size_y = int(math.floor(float(y_max - y_min + 1 + (n_reg_y - 1) * overlap_spat) / float(n_reg_y)));

        # Read template parameter file
        try:
            with open(template_file) as par_file:
                template_par = par_file.readlines()
        except:
            sys.stderr.write("Error: Failed to read template parameter file sofia.par.\n")
            sys.exit(16)

        sofia_params = []
        # Create set of parameter files
        for y in range(n_reg_y):
            for x in range(n_reg_x):
                index = x + (n_reg_x * y) + 1;

                par = template_par[:]

                x1 = x_min + int(math.floor(x * (size_x - overlap_spat)));
                y1 = y_min + int(math.floor(y * (size_y - overlap_spat)));
                x2 = x1 + size_x;
                y2 = y1 + size_y;

                if(x1 < 0): x1 = 0;
                if(y1 < 0): y1 = 0;
                if(x2 > x_max or x == n_reg_x - 1): x2 = x_max;
                if(y2 > y_max or y == n_reg_y - 1): y2 = y_max;

                # Add input region to par file
                i = self.substr_search(par, "input.region")
                if(i < 0): par.append("input.region  =  {0:d},{1:d},{2:d},{3:d},{4:d},{5:d}\n".format(x1, x2, y1, y2, z_min, z_max));
                else: par[i] = "input.region  =  {0:d},{1:d},{2:d},{3:d},{4:d},{5:d}\n".format(x1, x2, y1, y2, z_min, z_max);

                # Add output directory to par file
                i = self.substr_search(par, "output.directory")
                if(i < 0): par.append("output.directory  =  {0}\n".format(output_directory));
                else: par[i] = "output.directory  =  {0}\n".format(output_directory);

                # Add output filename to par file
                output_filename = '{0}_{1:03d}'.format(output_db_name, index)
                i = self.substr_search(par, "output.filename");
                if(i < 0): par.append("output.filename  =  {0}\n".format(output_filename));
                else: par[i] = "output.filename  =  {0}\n".format(output_filename);
    
                # Dump parameters into new file
                par_filename = "sofia_{0:03d}.par".format(index)
                sys.stdout.write("  {0}/{1}\n".format(par_file_directory, par_filename));
                try:
                    with open(par_filename, "w") as par_file:
                        for item in par:
                            par_file.write("{0}".format(item));
                except:
                    sys.stderr.write("Error: Failed to write output parameter file: {}\n".format(par_filename));
                    sys.exit(1);

                    sofia_params.append(SofiaParams(par_filename, par_file_directory, output_filename, output_directory))

        return sofia_params
