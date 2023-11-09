from pathlib import Path
import os
import time


def main():
    jobs_pth = Path("jobs")
    for job in jobs_pth.glob("*"):
        os.system("sh build.sh")
        os.system(f"bsub < {job}")
        time.sleep(10)

if __name__ == "__main__":
    main()