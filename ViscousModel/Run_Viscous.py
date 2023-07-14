import sys
import os

arg1 = sys.argv[1] # Cores
arg2 = sys.argv[2] # Test Directory
arg3 = sys.argv[3] # Test Name
testDir = f"{arg2}/Outputs"
subDir = testDir + "/"
testName = f"{arg2}/{arg3}"
os.system(f"rm -rf {testDir} && mkdir {testDir}")
os.system(f"mpirun -n {arg1} python3 {testName} {testDir} {subDir}")
