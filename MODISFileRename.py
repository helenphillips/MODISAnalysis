import os, sys

oldproduct = "_" + sys.argv[1] + ".asc"
newproduct = "___" + sys.argv[1] + ".asc"
directory = sys.argv[2]

for filename in os.listdir("."):
	if filename.endswith(product):
		x = filename.find(product)
		newname = filename[:x]
		newname = newname + newproduct
		os.rename(filename, newname)
