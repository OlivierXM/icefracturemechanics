#!/usr/bin/python3
"""
    Generate a packaged directory to copy to HPC
    All run files necessary for the test are in copyList
"""
import os
import sys

directoryName = input("Name directory:")
if (os.path.isfile(directoryName)):
    raise Exception("File by same name found!",f"Cannot create {directoryName} as a file by the same name exists!")
    
testName = input("Name TestFile:")
os.system(f"rm -rf {directoryName}")
os.system(f"mkdir {directoryName}")

# List of files to pack into the destination folder
copyList = ["Bar.h5","Bar.xdmf","../Utilities/Notes.py",
            "../Utilities/FileNamer.py", testName, "../Utilities/cdfx.py","ViscousMats.py",
            "Run_Viscous.py","./'Creep Testing'/LoadFile.dat","TestParameters.py"
           ]
        
for item in copyList:
    if (item == testName):
        os.system(f"cp {item} {directoryName}/ViscousTest.py")
    else:
        os.system(f"cp {item} {directoryName}/.")
## END SCRIPT ##
