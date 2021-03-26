# queue of processes
# parallel vs serial
# keep track of jobs, system resources.



#HPC cluster
# generate HPC job request
# trigger
# Keep track of process/job
from multiprocessing import Pool
import time

class SoFiA:
    PROCESS_NAME = 'sofia'
    @classmethod
    def execute(cls, config):
        time.sleep(2)
        with open(config, 'r') as conf_file:
            print(config)
    
class Executor:
    def __init__(self, max_parallel_process):
        self.max_parallel_process = max_parallel_process

    def run(self, configs):
        print("executing", configs)
        with Pool(processes=self.max_parallel_process) as pool:
            print(pool.map(SoFiA.execute, configs))
    
