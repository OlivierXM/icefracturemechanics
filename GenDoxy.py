#!/usr/bin/python3
"""

Use an existing dolfinx conda env to generate documentation

"""
import os

styleOut = "conda"

## Conda Details
condaName = "dolfinx_1"

## List of files to copy
copyList = ["ViscousModel", "Utilities",
            "Doxyfile", "README.md", "MarkDownDocuments"
           ]
os.system(f"rm -rf html && mkdir html")

if (styleOut == "conda"):
    os.system("conda init bash")
    os.system("rm -rf DoxyGen && mkdir DoxyGen")
    os.system(f"conda activate {condaName}")
    os.system("doxygen Doxyfile")
    os.system("cp -rf DoxyGen/html .")
    os.system("rm -rf DoxyGen")
    os.system(f"conda deactivate")


## END SCRIPT ##
