import os, sys

oldproduct = "_" + sys.argv[1] + ".asc"
newproduct = "___" + sys.argv[1] + ".asc"
directory = os.path.abspath(sys.argv[2])

for filename in os.listdir(directory):
	if filename.endswith(oldproduct):
		x = filename.find(oldproduct)
		newname = filename[:x]
		newname = newname + newproduct
		os.rename(os.path.join(directory,filename), os.path.join(directory, newname))
