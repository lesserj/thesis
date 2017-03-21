import os
import sys
import glob
import shutil
import csv
import ConfigParser
#import argparse

#input arguments should be in the form (path-row-year), (path-row-year)... output

#Load configuration information from settings.ini
config = ConfigParser.ConfigParser()
config.read("settings.ini")

training_data_path = config.get("Directories", "TRAINING_DATA_OUTPUT_DIRECTORY")
output_path = config.get("Directories", "COMBINED_TRAINING_DATA_OUTPUT_DIRECTORY") 
c5_exe = config.get("Executables", "C5_EXE") 
cubist_exe = config.get("Executables", "CUBIST_EXE") 

c5_names_template = config.get("Templates", "C5_NAMES_TEMPLATE")
cubist_names_template = config.get("Templates", "CUBIST_NAMES_TEMPLATE")

years = [2008,2009,2010]
scenes = []

#output_folder = ''
#for arg in sys.argv[1:-1]:
#	part = arg.split('-')
#	output_folder += '_' + part[0] + part[1] + part[2]

# We can really add a lot of options in a more robust manner here, but this will work for now.
if "-x8" in sys.argv:
	path = int(sys.argv[2])
	row = int(sys.argv[3])
	output_folder = str(path) + "_" + str(row)
	
	for year in years:
		scenes.append([path+1,row,year])
	for year in years:
		scenes.append([path+1,row+1,year])
	for year in years:
		scenes.append([path+1,row-1,year])
	for year in years:
		scenes.append([path,row+1,year])
	for year in years:
		scenes.append([path,row-1,year])
	for year in years:
		scenes.append([path-1,row,year])
	for year in years:
		scenes.append([path-1,row+1,year])
	for year in years:
		scenes.append([path-1,row-1,year])
	
else:
	raise SystemExit
	
output_path = os.path.join(output_path, output_folder)
os.mkdir(output_path)
shutil.copyfile(c5_names_template, os.path.join(output_path,'train_c5.names'))
shutil.copyfile(cubist_names_template, os.path.join(output_path,'train_cubist.names'))

c5_outfile = open(os.path.join(output_path,'train_c5.data'), 'w')

for scene in scenes:
	print scene
	inpath = os.path.join(training_data_path, str(scene[0]) + "_" + str(scene[1]), str(scene[2]))

	for file in glob.glob(inpath + '_*.txt'):
		with open(file) as infile:
			c5_outfile.write(infile.read())
print inpath			
c5_outfile.close()

with open(os.path.join(output_path,'train_c5.data'), "rb") as infile, open(os.path.join(output_path,'train_cubist.data'), "wb") as outfile:
	r = csv.reader(infile)
	
	for row in r:
		if row[-1] != '0':
			row[-1] = "100"
		outfile.write(",".join(row))
		outfile.write("\n")

os.system( cubist_exe + ' -f ' + os.path.join(output_path,'train_cubist') )		
#os.remove(os.path.join(output_path,'train_cubist.data'))
raise SystemExit	
os.system( c5_exe + ' -f ' + os.path.join(output_path,'train_c5') )

