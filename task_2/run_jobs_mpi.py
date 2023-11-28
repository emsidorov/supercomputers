from pathlib import Path
import os
import time
from math import pi
import math


def main():
    gs = [128, 256]
    ps = [1, 2, 4]
    ts = [1, 2, 4]
    Ls = [1.0]
    
    for grid_size in gs:
        for p in ps:
            for t in ts:
                for L in Ls:
                    for i in range(3):
                        file_out = f"results_mpi/{grid_size}_{p}_{t}_{int(L)}_{i}.out"
                        os.system(
                            f"mpisubmit.pl -p {p} -t {t} --stdout {file_out} --stderr tmp.err mpi {L} {L} {L} {grid_size} 0.4 20 0"
                        )
                        
                        time.sleep(60)
                                
                                


if __name__ == "__main__":
    main()