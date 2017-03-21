import os
import sys
import glob
import ConfigParser
import gdal
#import osr

#parameters
year = sys.argv[1]

#Load configuration information from settings.ini
config = ConfigParser.ConfigParser()
config.read("settings.ini")

baer_dir = config.get("Directories", "BAER_DATA_DIRECTORY")

data_dir = os.path.join(baer_dir,year)
print data_dir
perimeter_file = open(os.path.join(data_dir, 'perimeter.txt'), 'w')

# USGS Albers Projection
#PROJ = 'PROJCS["Albers_Conic_Equal_Area",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["meters",1]]'
#output_srs = osr.SpatialReference(PROJ)

#Find all files that are barc4
for fire in glob.glob(os.path.join(data_dir, '*', '*barc4_alb.img')):
	print fire
	
	#split file by the "_" separator and iterate through to identify date in file path
	file_path = fire.split('_')
	for date in file_path:
		if date.startswith(str(year)):
			barc4_ds = gdal.Open(fire, 0)
			
			#Calculate Extent
			geoTransform = barc4_ds.GetGeoTransform()
			x1 = geoTransform[0]
			x2 = geoTransform[0] + (barc4_ds.RasterXSize * geoTransform[1])
			y1 = geoTransform[3]
			y2 = geoTransform[3] + (barc4_ds.RasterYSize * geoTransform[5])
			
			#make sure all extents are calculated in the same Projection
			#(I've manually reprojected all the baer data to albers.  This
			# code is no longer needed, but left for historical purposes)
			#input_srs = osr.SpatialReference()
			#input_srs.ImportFromWkt(barc4_ds.GetProjectionRef())
			#transform = osr.CoordinateTransformation(input_srs,output_srs)
			#corner_nw = transform.TransformPoint(x1,y1)
			#corner_sw = transform.TransformPoint(x2,y2)
			corner_nw = (x1, y1)
			corner_sw = (x2, y2)

			#Output information to csv file as file path, extent, date
			perimeter_file.write(fire + "," + str(corner_nw[0]) + "," + str(corner_nw[1])  + "," + str(corner_sw[0]) + "," + str(corner_sw[1]) + "," + date + "\n")
			
			#We don't want to output the same file multiple times if the date appears in the file path multiple times
			break
perimeter_file.close()