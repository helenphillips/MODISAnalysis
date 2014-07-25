MODISFileRename.py
===================


Python code to rename old MODIS subsets file names into the new format.

Example:
Old file name: lat_long_enddate_MOD13Q1.asc
New file name: lat_long_enddate___MOD13Q1.asc

How to use:
-----------
- Save the file to a directory on your computer.
- Open terminal/command line, and navigate to the .py file.
- Type 'python MODISFileRename.py [product] [folder]

product: string; the product name given in the current file name (example: MOD13Q1, or if two products, MCD43A4_MCD43A2).

folder: string; path to the folder where subsets are downloaded ("." if current directory)


This has been tested using python 2.7, please let me know if it doesn't work with other versions.