import os
import sys
import glob
import math
import random
import ConfigParser
import subprocess
from datetime import date
import random
import gdal
import ogr
import osr
#import numpy as np

# USGS Albers Projection
PROJ = 'PROJCS["Albers_Conic_Equal_Area",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["meters",1]]'

# Command line arguments to determine path, row, and year of data to extract
path = int(sys.argv[1])
row = int(sys.argv[2])
year = int(sys.argv[3])

#Parse configuration information from the settings.ini file
config = ConfigParser.ConfigParser()
config.read("settings.ini")

#Load configuration information from settings.ini
landsat_data_directory  = config.get("Directories", "LANDSAT_DATA_DIRECTORY")
baer_dir = config.get("Directories", "BAER_DATA_DIRECTORY")
training_data_output_directory = config.get("Directories", "TRAINING_DATA_OUTPUT_DIRECTORY")
fire_history_path = config.get("Directories", "FIRE_HISTORY_DIRECTORY")

# Processed Landsat data location
data_dir = os.path.join(landsat_data_directory, "p0" + str(path) + "r0" + str(row), 'HFA')
baer_dir = os.path.join(baer_dir,str(year))
outfile_path = os.path.join(training_data_output_directory, str(path) + "_" + str(row), str(year))

try:
	os.makedirs(os.path.dirname(outfile_path) + "\\")
except WindowsError:
	#directory already exists or other permissions/disk errors
	pass
except OSError:
	#Other operating system, and windows write errors
	pass

# This is date order, but not necessarily leaf-on/leaf-off order.
# Get three years of data (two scenes per year).  Middle scene represents
# the year from the command line argument and should represent the 
# year of disturbance
years = [year-1,year,year+1]

# Only extract data where disturbance raster is one of these values
# Baer data has the following classes.  We are only likely to pick up class 3 and 4
#	0 = outside perimeter
#	1 = unchanged / very low (Dark Green) | This means the area after the fire was indistinguishable from pre-fire conditions. This does not always indicate the area did not burn.
#	2 = low severity (Cyan) | This severity class represents areas of surface fire with little change in cover and little mortality of the dominant vegetation.
#	3 = moderate severity (Yellow) | This severity class is between low and high and means there is a mixture of effects on the dominant vegetation.
#	4 = high severity (Red) | This severity class represents areas where the dominant vegetation has high to complete mortality.
disturbance_values = [2,3,4]

#Create separate text files for each disturbance type and for undisturbed data
outfiles = []
for disturbance_value in disturbance_values:
	outfiles.append(open(outfile_path + "_" + str(disturbance_value) + '.txt', 'w'))

#undisturbed values
outfiles.append(open(outfile_path + "_0" + '.txt', 'w')) 

# Tracking coordinates this way and extracting pixels not in the coordinates set is slow as molasses
# Try using a numpy mask and drawing points randomly could drastically speed up this operation
coordinates = set() #coordinates of all baer pixels

# Lists to store rasters and raster metadata
rasters = []
scene_dates = []
days_of_year = []

# Identify Landsat rasters for analysis
for raster_path in glob.glob(os.path.join(data_dir, '*.img.gz')):
	#Only get the Landsat scene rasters.  We could do something more elegant with regular
	#expressions, but this gets the job done well enough for now
	if len(os.path.basename(raster_path)) == 27:
		#only get rasters within the year range we are extracting from
		if int(os.path.basename(raster_path)[12:16]) in years:
			#Open rasters while still gzipped.  This is slower to read, but probably faster and
			#cleaner than extracting first
			#rasters.append(gdal.Open('/vsigzip/' + raster_path ,0))
			if not os.path.isfile(raster_path[0:-3]):
				proc = subprocess.Popen(['7-Zip\\7z.exe', '-o' + os.path.dirname(raster_path), 'e', raster_path], stdout=subprocess.PIPE, shell=True)
				out, err = proc.communicate()
				if err:
					print raster_path[0:-3]
					raise SystemExit
				proc, out, err = None, None, None
			rasters.append(gdal.Open(raster_path[0:-3] ,0))
			#Extract date of scene from raster file path and store in a list
			scene_dates.append(int(os.path.basename(raster_path)[12:20]))
	raster_path = None
	
if len(rasters) > 10:
	print ("ERROR: More than 10 Rasters")
	
	for raster_path in glob.glob(os.path.join(data_dir, '*.img.gz')):
		if os.path.isfile(raster_path[0:-3]): 
			os.remove(raster_path[0:-3])
			
	raise SystemExit
	
# Calculate the day of year each scene was taken using the scene dates
for scene_date in scene_dates:
	scene_year = int(str(scene_date)[0:4])
	scene_month = int(str(scene_date)[4:6])
	scene_day = int(str(scene_date)[6:8])
	scene_doy = date(scene_year,scene_month,scene_day).timetuple().tm_yday
	days_of_year.append(scene_doy)

# For now we'll assume all rasters in a scene have the same extents/projections
# and extract the geotransform data from the first raster
# At some point if that changes, then we will have to project into a uniform projection
# (if different) and compute the extent as the area where all scenes have data
scene_geotransform = rasters[0].GetGeoTransform() 
scene_upper_left_x = scene_geotransform[0]
scene_pixel_width = scene_geotransform[1]
scene_upper_left_y = scene_geotransform[3]
scene_pixel_height = scene_geotransform[5]
scene_xsize = rasters[0].RasterXSize
scene_ysize = rasters[0].RasterYSize

# loads a list of our rasters extent in the order [left,top,right,bottom]
extent = [
	scene_upper_left_x, 
	scene_upper_left_y, 
	scene_upper_left_x + (scene_xsize * scene_pixel_width), 
	scene_upper_left_y + (scene_ysize * scene_pixel_height)
	]

# Create polygon of scene extent to intersect with baer/disturbance
# data to identify disturbances overlapping the scenes
scene_geometry = ogr.Geometry(ogr.wkbPolygon)
ring = ogr.Geometry(ogr.wkbLinearRing)
ring.AddPoint(extent[0],extent[3]) #x1y1
ring.AddPoint(extent[2],extent[3]) #x2y2
ring.AddPoint(extent[2],extent[1]) #x2y1
ring.AddPoint(extent[0],extent[1]) #x1y2
ring.CloseRings()
scene_geometry.AddGeometry(ring)
ring = None

shp_driver = ogr.GetDriverByName("ESRI Shapefile")
fire_history_ds = shp_driver.Open(os.path.join(fire_history_path, str(year) + '.shp'))
print os.path.join(fire_history_path, str(year) + '.shp')
if fire_history_ds is None:
	print 'could not open fire history perimeters %s' % os.path.join(fire_history_path, str(year) + '.shp')
	raise SystemExit
	
fire_geometries = []	
fire_history_lyr = fire_history_ds.GetLayer()
fire_history_lyr.SetSpatialFilter(scene_geometry)


# numpy array to mask out disturbed pixels which we'll use to then extract a
# sample of undisturbed pixels
#mask = np.ones((scene_xsize,scene_ysize))

# Keeps track of the number of disturbed pixels
count = 0

#Read in fire extents and see if they intersect the scene geometry
with open(os.path.join(baer_dir,'perimeter.txt')) as f:
	for line in f:
		# last character is a newline character, so we trim that off
		# and split into a comma separated list
		baer =  line[:-1].split(',')
		baer_file_path = baer[0]
		baer_left = float(baer[1])
		baer_top = float(baer[2])
		baer_right = float(baer[3])
		baer_bottom = float(baer[4])
		baer_date = int(baer[5]) #yyyymmdd
		baer = None
		
		#Check if fire occurs between 3rd and 5th Landsat scenes in the stack
			#  ____	  ____   ____   ____   ____	  ____
			# |	y1 | | y1 | | y2 | | y2 | | y3 | | y3 |
			# | l1 | | l2 | | l1 | | l2 | | L1 | | l2 |
			# |____| |____| |____| |____| |____| |____|
			#                     ^      ^
			#					  | fire | 
		######################################################################
		
		if baer_date > scene_dates[2] and baer_date < scene_dates[4]:
			data = []
				
			# Create polygon of the baer/disturbance extent
			baer_geometry = ogr.Geometry(ogr.wkbPolygon)
			ring = ogr.Geometry(ogr.wkbLinearRing)
			ring.AddPoint(baer_left,baer_bottom) #x1y1
			ring.AddPoint(baer_right,baer_bottom) #x2y1
			ring.AddPoint(baer_right,baer_top) #x2y2
			ring.AddPoint(baer_left,baer_top) #x1y2
			ring.CloseRings()
			baer_geometry.AddGeometry(ring)
			ring = None

			# If the baer/disturbance geometry intersects the scene geometry
			# then the scene likely has been disturbed so extract disturbed pixels
			if scene_geometry.Intersects(baer_geometry):
				# Open baer/disturbance raster as read only
				baer_ds = gdal.Open(baer_file_path,0)
				baer_geotransform = baer_ds.GetGeoTransform()
				
				baer_upper_left_x = baer_geotransform[0]
				baer_pixel_width = baer_geotransform[1]
				baer_we_rotation = baer_geotransform[2]
				baer_upper_left_y = baer_geotransform[3]
				baer_ns_rotation = baer_geotransform[4]
				baer_pixel_height = baer_geotransform[5]
				baer_xsize = baer_ds.RasterXSize
				baer_ysize = baer_ds.RasterYSize
				
				# Sometimes BAER data is in UTM projection instead of Albers projection
				# Project BAER data into Albers, if necessary, using Nearest Neighbor algorithm
				print baer_ds.GetProjection()
				raise SystemExit
				if "utm" in baer_file_path:
					print 'utm'
					# Define transformation from input projection to output projection
					utm_proj = osr.SpatialReference(baer_ds.GetProjection())
					alb_proj = osr.SpatialReference(PROJ)
					transformation = osr.CoordinateTransformation(utm_proj, alb_proj)
					
					# Calculate geotransform of projected raster
					# top left corner
					(ulx, uly, ulz) = transformation.TransformPoint(baer_upper_left_x, baer_upper_left_y)
					# bottom right corner
					(lrx, lry, lrz) = transformation.TransformPoint(baer_upper_left_x + baer_pixel_width*baer_xsize, baer_upper_left_y + baer_pixel_height*baer_ysize)
					# bottom left corner
					(llx, lly, llz) = transformation.TransformPoint(baer_upper_left_x, baer_upper_left_y + baer_pixel_height*baer_ysize)
					# top right corner
					(urx, ury, urz) = transformation.TransformPoint(baer_upper_left_x + baer_pixel_width*baer_xsize, baer_upper_left_y)

					transformation = None
					
					# Identify the furthest left point and the top most point
					# Sometimes the rasters in a particular projection are skewed relative to the new position
					# in a way that makes the bottom left corner "more left" than the top corner, and similar
					# with the top and bottom most values.  By detecting the left most and top most corners,
					# we ensure that we output a new projected raster with the correct dimensions.
					lx = min(llx, ulx)
					rx = max(urx, lrx)
					ty = max(uly, lly)
					by = min(lry, ury)
					
					# Create empty raster in new projection with the correct bounds
					mem = gdal.GetDriverByName('MEM')
					dest = mem.Create('', int((rx - lx)/baer_pixel_width), int((by - ty)/ baer_pixel_height), 1, baer_ds.GetRasterBand(1).DataType)
					new_geotransform = ( int(lx), baer_pixel_width, baer_we_rotation, int(ty), baer_ns_rotation, baer_pixel_height )
					dest.SetGeoTransform( new_geotransform )
					new_geotransform = None
					lr,rx,ty,by = None,None,None,None
					dest.SetProjection ( alb_proj.ExportToWkt() )
					
					# Reproject the image from the original dataset into the new dataset using Nearest Neighbor
					# We use nearest neighbor because it is most important to preserve actual pixel values
					# rather than altering them with an algorithm that "averages" values
					gdal.ReprojectImage( baer_ds, dest, utm_proj.ExportToWkt(), alb_proj.ExportToWkt(), gdal.GRA_NearestNeighbour  )
					utm_proj, alb_proj = None, None
					
					# Save raster to disk
					#driver = gdal.GetDriverByName ( "GTiff" )
					#dst_ds = driver.CreateCopy( os.path.join("c:/",'temp/', baer[5] + r'.tif'), dest, 0 )
					#dst_ds = None # Flush the dataset to disk
					
					baer_ds = dest
					baer_geotransform = baer_ds.GetGeoTransform()
					baer_upper_left_x = baer_geotransform[0]
					baer_pixel_width = baer_geotransform[1]
					baer_we_rotation = baer_geotransform[2]
					baer_upper_left_y = baer_geotransform[3]
					baer_ns_rotation = baer_geotransform[4]
					baer_pixel_height = baer_geotransform[5]
					baer_xsize = baer_ds.RasterXSize
					baer_ysize = baer_ds.RasterYSize
					
					# The geometry of the new reprojected raster extent is different and so needs to be recalculated
					baer_geometry = None
					baer_geometry = ogr.Geometry(ogr.wkbPolygon)
					ring = ogr.Geometry(ogr.wkbLinearRing)
					ring.AddPoint(ulx, uly) #x1y1
					ring.AddPoint(urx, ury) #x2y1
					ring.AddPoint(lrx, lry) #x2y2
					ring.AddPoint(llx, lly) #x1y2
					ring.CloseRings()
					baer_geometry.AddGeometry(ring)
					ring, ulx, uly, urx, ury, lrx, lry, llx, lly = None, None, None, None, None, None, None, None, None
					
				disturbance_extent = scene_geometry.Intersection(baer_geometry).GetEnvelope()
				disturbance_min_x = disturbance_extent[0]
				disturbance_max_x = disturbance_extent[1]
				disturbance_min_y = disturbance_extent[2]
				disturbance_max_y = disturbance_extent[3]

				
				# Values to extract from Scene Data
				# NOTE: BAER imagery is not always snapped to LANDSAT imagery.  This can cause a 1 pixel
				# shift between baer values and underlying Landsat values
				# Since BAER data is for fires over 40 acres, a few edge pixels shouldn't make too big of 
				# a difference (sensitivity analysis)?  This shouldn't cause a problem detecting smaller 
				# fires.
				# on a side-note: do smaller fires have different spectral properties from big fires?
				if (abs((scene_upper_left_x - baer_upper_left_x) % scene_pixel_width) < (scene_pixel_width / 2)):
					x_snap_factor = 1
				else:
					x_snap_factor = 0
				
				if (abs((scene_upper_left_y  - baer_upper_left_y) % scene_pixel_height) < (scene_pixel_height/2)):
					y_snap_factor = -1
				else :
					y_snap_factor = 0

				print "snap factors: ", x_snap_factor,y_snap_factor
			
				# Subtract 1 from both x_size and y_size since we are rounding up the starting locations
				# (x_start, y_start) using the math.ceil function.  This ensures that we don't run out
				# of bounds on the East side of our Landsat rasters
				ls_x_start = int((disturbance_min_x - scene_upper_left_x) / scene_pixel_width) + x_snap_factor
				ls_x_size = int((disturbance_max_x - disturbance_min_x) / scene_pixel_width) - x_snap_factor
				ls_y_start = int((disturbance_max_y - scene_upper_left_y) / scene_pixel_height) + y_snap_factor
				ls_y_size = int((disturbance_min_y - disturbance_max_y) / scene_pixel_height) - y_snap_factor
				
				for raster in rasters:
					for band in range(raster.RasterCount):
						data.append(raster.GetRasterBand(band+1).ReadAsArray(ls_x_start, ls_y_start, ls_x_size, ls_y_size))
				
				#Values to extract from BAER
				baer_x_start = int((disturbance_min_x - baer_upper_left_x) / baer_pixel_width)
				baer_x_size = int((disturbance_max_x - disturbance_min_x) / baer_pixel_width)
				baer_y_start = int((disturbance_max_y - baer_upper_left_y) / baer_pixel_height)
				baer_y_size = int((disturbance_min_y - disturbance_max_y) / baer_pixel_height)
				
				if baer_x_start < 0:
					baer_x_start = 0
					baer_x_size -= 1
					
				if baer_y_start < 0:
					baer_y_start = 0
					baer_y_size -= 1
					
				data.append(baer_ds.GetRasterBand(1).ReadAsArray(baer_x_start, baer_y_start, baer_x_size, baer_y_size))
				
				for column in xrange(min(ls_x_size,baer_x_size)-1):
					for row in xrange(min(ls_y_size,baer_y_size)-1):
						l=[]
						
						x = str(scene_upper_left_x  + (0.5 * scene_pixel_width) + ((ls_x_start + column) * scene_pixel_width))
						y = str(scene_upper_left_y + (0.5 * scene_pixel_height) + ((ls_y_start + row) * scene_pixel_height))
						l.append(x)
						l.append(y)
						
						coordinates.add((int(float(x)),int(float(y))))
						
						for raster in xrange(len(data) - 1):
							l.append(str(data[raster][row][column]))
							
						for doy in days_of_year:
							l.append(str(doy))
							
						l.append(str(data[len(data) - 1][row][column]))
						
						# If any pixel has no value in the Landsat stack, then ignore the entire pixel stack
						if all(i != '0' for i in l[2:-6]) and (int(l[-1:][0]) in disturbance_values):
							outfiles[disturbance_values.index(int(l[-1:][0]))].write((',').join(l) + '\n')
							count += 1
							
						l = None
		
				baer_x_start = None
				baer_x_size = None
				baer_y_start = None
				baer_y_size = None
		
			data = None
		
		baer_ds = None 
		baer_date = None
		baer_file_path = None
		baer_left = None
		baer_top = None
		baer_right = None
		baer_bottom = None
		baer_date = None
	
	line = None
		
# read undisturbed points

max_undisturbed = min(count*10, 1000000)
x_size = rasters[0].RasterXSize
y_size = rasters[0].RasterYSize

pixels = x_size * y_size
try:
	offset = pixels / max_undisturbed
except ZeroDivisionError:
	print 'no pixels to extract'
	for outfile in outfiles:
		outfile.close()
	
	rasters = None
	
	for raster_path in glob.glob(os.path.join(data_dir, '*.img.gz')):
		if os.path.isfile(raster_path[0:-3]): 
			pass
			#os.remove(raster_path[0:-3])
	
	raise SystemExit #Don't process undisturbed scenes with no disturbance pixels
	#offset = 1000
print pixels
print offset

#read one line at a time, take pixel sample by offset
for row in random.sample(range(0,y_size), int(math.sqrt(max_undisturbed))):

	data = []
	for raster in rasters:
		for band in range(raster.RasterCount):
			data.append(raster.GetRasterBand(band+1).ReadAsArray(0, row, x_size, 1))

	for column in random.sample(range(0,x_size), int(math.sqrt(max_undisturbed))):
				
		x = str(scene_geotransform[0] + (0.5 * scene_geotransform[1]) + (column * scene_geotransform[1]))
		y = str(scene_geotransform[3] + (0.5 * scene_geotransform[5]) + (row * scene_geotransform[5]))
		if (int(float(x)),int(float(y))) in coordinates: continue
		
		in_perimeter = False
		fire_history_lyr.ResetReading()
		for fire_perimeter in fire_history_lyr:
			geom = fire_perimeter.GetGeometryRef()
			pt = ogr.Geometry(ogr.wkbPoint)
			pt.SetPoint(0, float(x), float(y))
			# test
			if pt.Within(geom):
				#print fire_perimeter.GetFieldAsString(7), x, y
				in_perimeter = True
				break
		
		if in_perimeter == True: 
			continue

		l=[]
		l.append(x)
		l.append(y)
		
		for raster in xrange(len(data)):
			l.append(str(data[raster][0][column]))

		for doy in days_of_year:
			l.append(str(doy))	
			
		l.append('0') #undisturbed class
		if all(i != '0' for i in l[2:-7]):
			outfiles[-1].write((',').join(l) + '\n')

for outfile in outfiles:
	outfile.close()
	
rasters = None
	
for raster_path in glob.glob(os.path.join(data_dir, '*.img.gz')):
	if os.path.isfile(raster_path[0:-3]): 
		pass
		#os.remove(raster_path[0:-3])