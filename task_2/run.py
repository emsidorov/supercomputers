import os
from math import pi
import math


def main():
    gs = [128, 256]
    ts = [1, 2, 4, 8, 16, 32]
    Ls = [1.0, pi]
    
    for grid_size in gs:
        for t in ts:
            for L in Ls:
                L_name = int(L)
                p = int(math.ceil(t / 8.0))
                job_text = f'#BSUB -n {p}\n#BSUB -x\n#BSUB -W 00:15\n#BSUB -o "results/{grid_size}_{L_name}_{t}.out"\n#BSUB -e "tmp.err"\n#BSUB -R "span[hosts=1]"\nOMP_NUM_THREADS={t} ./omp {L} {L} {L} {grid_size} 0.01 20'
                job_name = f"jobs/{grid_size}_{L_name}_{t}.lsf"
                with open(job_name, "w") as file:
                    file.write(job_text)

if __name__ == "__main__":
    main()