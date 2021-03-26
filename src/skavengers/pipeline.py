# Create Configs
# Launch Process
# combine result

import sys
import configparser

# config class

class Pipeline:
    def __init__(self, pipeline_setup, config_file):
        self.config = configparser.ConfigParser()
        success = self.config.read(config_file)

        if(len(success) == 0):
            sys.stderr.write("Error: Failed to read config file: pipeline_setup.ini\n")
            sys.exit(1)

        self.pipeline_setup = pipeline_setup
        self.config_path = config_path
        self.sofia_params = []

    def init(self):
        self.sofia_params = self.pipeline_setup.gen_par_files(self.config)

    def execute(self, executor):
        executor.run([sofia_param.par_file for sofia_param in self.sofia_params])


            
    def finish(self):
        #merge catalogue results
        pass


