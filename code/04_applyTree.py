#applyTree
#	applytree.py [path] [row] [year]
#	apply tree to set of scenes

import os
import sys
import glob
import shutil
import struct
import subprocess
import ConfigParser
import numpy as np
from datetime import date
try:
	from osgeo import gdal
except ImportError:
	import gdal

#http://stackoverflow.com/questions/4914008/efficient-way-of-parsing-fixed-width-files-in-python
def slices(s, *args):
    position = 0
    for length in args:
        yield s[position:position + length]
        position += length
		
path = int(sys.argv[1])
row = int(sys.argv[2])
year = int(sys.argv[3])

#Load configuration information from settings.ini
config = ConfigParser.ConfigParser()
config.read("settings.ini")

landsat_data_path = config.get("Directories", "LANDSAT_DATA_DIRECTORY")
combined_training_data_output_directory = config.get("Directories", "COMBINED_TRAINING_DATA_OUTPUT_DIRECTORY")
cubist_interpreter_exe = config.get("Executables", "CUBIST_INTERPRETER_EXE")

# for each row in set of rasters
	#export row data to text file of cases
	#process cases with c5
	#import results into new raster
	
# Processed Landsat data location
data_dir = os.path.join(landsat_data_path, "p0" + str(path) + "r0" + str(row), 'HFA')

# Output location
output_dir = os.path.join(combined_training_data_output_directory, str(path) + "_" + str(row))
output_ds_path = os.path.join(combined_training_data_output_directory, str(path) + "_" + str(row), str(year) + ".img")


try:
	os.makedirs(output_dir)
except WindowsError:
	#directory already exists
	pass

# This is date order, but not necessarily leaf-on/leaf-off order.
# Get three years of data (two scenes per year).  Middle scene represents
# the year from the command line argument and should represent the 
# year of disturbance.  This needs to match the years extracted as part
# of extract.py
years = [year-1,year,year+1]

scene_dates = []
days_of_year = []

# Identify Landsat rasters for analysis
rasters = []
for raster_path in glob.glob(os.path.join(data_dir, '*.img.gz')):
	#Only get the Landsat scene rasters.  We could do something more elegant with regular
	#expressions, but this gets the job done well enough for now
	if len(os.path.basename(raster_path)) == 27:
		#only get rasters within the year range we are extracting from
		if int(os.path.basename(raster_path)[12:16]) in years:
			#Open rasters while still gzipped.  This is slower to read, but probably faster and
			#cleaner than extracting first
			if not os.path.isfile(raster_path[0:-3]):
				proc = subprocess.Popen(['7-Zip\\7z.exe', '-o' + os.path.dirname(raster_path), 'e', raster_path], stdout=subprocess.PIPE, shell=True)
				out, err = proc.communicate()
				if err:
					print raster_path[0:-3]
					raise SystemExit
				proc, out, err = None, None, None
			rasters.append(gdal.Open(raster_path[0:-3] ,0))
			scene_dates.append(int(os.path.basename(raster_path)[12:20]))

if len(rasters) > 10:
	print ("ERROR: More than 10 Rasters")
	raise SystemExit
			
for scene_date in scene_dates:
	year = int(str(scene_date)[0:4])
	month = int(str(scene_date)[4:6])
	day = int(str(scene_date)[6:8])
	doy = date(year,month,day).timetuple().tm_yday
	days_of_year.append(doy)
	scene_date,year,month,day,doy = None,None,None,None,None

	
#Get raster properties
#This assumes all rasters have identical extents and projections
driver = rasters[0].GetDriver()
geoTransform = rasters[0].GetGeoTransform() 
projection = rasters[0].GetProjection()
xsize = rasters[0].RasterXSize
ysize = rasters[0].RasterYSize		
			
#Open output raster
output_ds = driver.Create(output_ds_path, xsize, ysize, 1, gdal.GDT_Byte)
output_ds.SetGeoTransform(geoTransform)
output_ds.SetProjection(projection)
output_band = output_ds.GetRasterBand(1)

#Process raster data
filestem = os.path.join(output_dir,'train_cubist')

#Create day of year string outside of loop, since it should be constant
doys = []
for doy in days_of_year:
	doys.append(str(doy))
doys_str = (',').join(doys)

#Get number of datapoints
datapoints = 0;
for raster in rasters:
	datapoints += raster.RasterCount

#The first two datapoints are the coordinates
datapoints += 2

#The second to last set of datapoints are place holders for the year
datapoints += len(days_of_year)
	
#The last datapoint is our unknown disturbance value labeled "?"	
datapoints += 1



for row in xrange(0,ysize,100):
	sys.stdout.write("%d\\%d  \r" % (row,ysize))
	
	if ((row + 100) > ysize):
		num_rows = ysize - row
	else:
		num_rows = 100

	data = np.zeros((num_rows * xsize * datapoints), dtype="S10")
	coordsX = np.zeros((num_rows * xsize), dtype="S10")
	coordsY = np.zeros((num_rows * xsize), dtype="S10")
	
	# Generate the coordinate values.
	#TODO: Put directly into our output array instead of an intermediate array
	i = 0
	for r in xrange(num_rows):
		for column in xrange(xsize):
			coordsX[i] = str(geoTransform[0] + (0.5 * geoTransform[1]) + (column * geoTransform[1]))
			coordsY[i] = str(geoTransform[3] + (0.5 * geoTransform[5]) + ((row + r) * geoTransform[5]))
			i += 1
	
	data[0::datapoints] = coordsX
	data[1::datapoints] = coordsY
	
	#Start at position 2 because pos 0 and 1 are the coords
	position = 2
	for raster in rasters:
		for band in range(raster.RasterCount):
			data[position::datapoints] = raster.GetRasterBand(band+1).ReadAsArray(0, row, xsize, num_rows).flatten()
			position += 1
	
	for doy in days_of_year:
		data[position::datapoints] = np.array([doy])
		position += 1
			
	data[position::datapoints] = np.array(["?"])

	np.savetxt(filestem + '.cases', data.reshape((num_rows * xsize, datapoints)), "%s", delimiter=",")
		
	proc = subprocess.Popen([cubist_interpreter_exe, '-f', filestem], stdout=subprocess.PIPE, shell=True)
	out, err = proc.communicate()

	outline = ''
	count = 0
	lines = out.split('\n')
	for line in lines:
		if line !="predicted values:" and line !="":
			outline += struct.pack('B',int(float(line)))
		count +=1

	output_band.WriteRaster( 0, row, xsize, num_rows, outline )   

rasters = None	
for raster_path in glob.glob(os.path.join(data_dir, '*.img.gz')):
	if os.path.isfile(raster_path[0:-3]): 
		#os.remove(raster_path[0:-3])
		pass

output_band.FlushCache()
output_band = None
output_ds = None