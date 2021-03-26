# queue of processes
# parallel vs serial
# keep track of jobs, system resources.



#HPC cluster
# generate HPC job request
# trigger
# Keep track of process/job
from multiprocessing import Pool
import time
import subprocess

class SoFiA:
    PROCESS_NAME = 'sofia'
    @classmethod
    def execute(cls, config):
        command = [cls.PROCESS_NAME, config]
        return subprocess.run(command, capture_output=True)
        
class Executor:
    def __init__(self, max_parallel_process):
        self.max_parallel_process = max_parallel_process

    def run(self, configs):
        completed_process = []
        with Pool(processes=self.max_parallel_process) as pool:
            completed_process = pool.map(SoFiA.execute, configs)

        return completed_process
