import os
import sys
import glob
import ConfigParser
import subprocess
import ogr
import gdal
import numpy as np

# Get command line variables
path = int(sys.argv[1])
row = int(sys.argv[2])
year = int(sys.argv[3])

#Load configuration information from settings.ini
config = ConfigParser.ConfigParser()
config.read("settings.ini")

combined_training_data_output_directory = config.get("Directories", 
                                "COMBINED_TRAINING_DATA_OUTPUT_DIRECTORY")
landsat_data_directory  = config.get("Directories", "LANDSAT_DATA_DIRECTORY")
data_dir = os.path.join(landsat_data_directory, 
                        "p0" + str(path) + "r0" + str(row), 'HFA')
    
###############################################################################
# Create data mask where all input Landsat pixels are not NoData
###############################################################################


first_ls_ds = gdal.Open('/vsigzip/' + glob.glob(os.path.join(data_dir, 
    '*.img.gz'))[0], 0)
driver = first_ls_ds.GetDriver()
landsat_xsize = first_ls_ds.RasterXSize
landsat_ysize = first_ls_ds.RasterYSize
geoTransform = first_ls_ds.GetGeoTransform()
proj = first_ls_ds.GetProjection()

data_mask_data = np.ones([landsat_ysize, landsat_xsize])

for landsat_path in glob.glob(os.path.join(data_dir, '*.img.gz')):
    print landsat_path
    landsat_ds = gdal.Open('/vsigzip/' + landsat_path, 0)
    #Assume NoData if band 1 is 0.  This is probably always, mostly true
    landsat_data = landsat_ds.GetRasterBand(1).ReadAsArray(0,0)
    landsat_data[landsat_data > 0] = 1
    data_mask_data = data_mask_data * landsat_data
    
data_mask_path = os.path.join(combined_training_data_output_directory, 
                                str(path) + "_" + str(row), 
                                str(year) + "_data_mask.img")

data_mask_ds = driver.Create(data_mask_path, 
                          landsat_xsize, landsat_ysize, 1, 1)
data_mask_ds.SetGeoTransform(geoTransform)
data_mask_ds.SetProjection(proj)
data_mask_ds.GetRasterBand(1).WriteArray(data_mask_data)


# This is not a pixel perfect vector as a lot of interpolation happens, 
# but it is good enough for our puproses of selecting fire points
output_vect_path = os.path.join(combined_training_data_output_directory, 
                str(path) + "_" + str(row), str(year) + "_data_mask.shp")
proc = subprocess.Popen(['C:\Program Files\GDAL\gdal_polygonize.py', 
                         '-8', 
                         '-f', 'ESRI Shapefile',
                         '-mask',data_mask_path, 
                         data_mask_path, 
                         output_vect_path], 
                         stdout=subprocess.PIPE, shell=True)
out, err = proc.communicate()
        
###############################################################################
# Reclassify map cubist output with everything below the detection threshold
# as 0 and everything equal to or above the detection threshold as 1
###############################################################################

DETECTION_THRESHOLD = 90

detected_disturbance_path = os.path.join(
    combined_training_data_output_directory, 
    str(path) + "_" + str(row), str(year) + ".img")

detection_mask_path = os.path.join(
    combined_training_data_output_directory, 
    str(path) + "_" + str(row), str(year) + "_detection_mask.img")

# Read in Raster
detected_disturbance_ds = gdal.Open(detected_disturbance_path,0)

# Get in raster properties
driver = detected_disturbance_ds.GetDriver()
detected_disturbance_band = detected_disturbance_ds.GetRasterBand(1)
detected_disturbance_data = detected_disturbance_band.ReadAsArray(0,0)

# reclass based on detection threshold
detected_disturbance_data[detected_disturbance_data 
    < DETECTION_THRESHOLD] = 0
detected_disturbance_data[detected_disturbance_data 
    >= DETECTION_THRESHOLD] = 1
    
# Mask NoData
detected_disturbance_data[data_mask_data == 0] = 255
data_mask_ds.GetRasterBand(1).SetNoDataValue(255)

# Write out new raster
detection_mask_ds = driver.CreateCopy(detection_mask_path, 
                                      detected_disturbance_ds, 0)
detection_mask_ds.GetRasterBand(1).WriteArray(detected_disturbance_data)
detection_mask_ds = None

# Sieve data to remove noise
destination_path = destination_path = os.path.join(
    combined_training_data_output_directory, 
    str(path) + "_" + str(row), str(year) + "_detection_mask-x8-11.img")
     
#11 pixels ~ 1 hectare     
proc = subprocess.Popen(['C:\Program Files\GDAL\gdal_sieve.py', 
                         '-st', '11', 
                         '-8',
                         detection_mask_path,'-of', 'HFA', destination_path
                         ], 
                         stdout=subprocess.PIPE, shell=True)
                             
out, err = proc.communicate()


###############################################################################
# Mosaic BAER data in path/row/year where the discover date is less than the
# last detection year scene date
###############################################################################

baer_dir = config.get("Directories", "BAER_DATA_DIRECTORY")
baer_dir = os.path.join(baer_dir,str(year))


scene_path = glob.glob(os.path.join(data_dir, 
                                    '*' + str(year) + '????.img.gz'))[1]
scene_date = scene_path.split('_')[2][7:11]

scene_ds = gdal.Open('/vsigzip/' + scene_path, 0)
scene_geotransform = scene_ds.GetGeoTransform() 
scene_upper_left_x = scene_geotransform[0]
scene_pixel_width = scene_geotransform[1]
scene_upper_left_y = scene_geotransform[3]
scene_pixel_height = scene_geotransform[5]
scene_xsize = scene_ds.RasterXSize
scene_ysize = scene_ds.RasterYSize

# loads a list of our rasters extent in the order [left,top,right,bottom]
extent = [
    scene_upper_left_x, 
    scene_upper_left_y, 
    scene_upper_left_x + (scene_xsize * scene_pixel_width), 
    scene_upper_left_y + (scene_ysize * scene_pixel_height)
    ]
    
scene_geometry = ogr.Geometry(ogr.wkbPolygon)
ring = ogr.Geometry(ogr.wkbLinearRing)
ring.AddPoint(extent[0],extent[3]) #x1y1
ring.AddPoint(extent[2],extent[3]) #x2y2
ring.AddPoint(extent[2],extent[1]) #x2y1
ring.AddPoint(extent[0],extent[1]) #x1y2
ring.CloseRings()
scene_geometry.AddGeometry(ring)
ring = None

baer_paths = []
#Read in fire extents and see if they intersect the scene geometry
with open(os.path.join(baer_dir,'perimeter.txt'), 'r') as f:
    for line in f:
        baer =  line[:-1].split(',')
        baer_file_path = baer[0]
        baer_left = float(baer[1])
        baer_top = float(baer[2])
        baer_right = float(baer[3])
        baer_bottom = float(baer[4])
        baer_date = baer[5][4:] #mmdd
        baer = None
        
        if baer_date <= scene_date:
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

            # If the baer/disturbance geometry intersects the scene 
            #geometry then the scene likely has been disturbed so 
            #extract disturbed pixels
            if scene_geometry.Intersects(baer_geometry):
                print baer_file_path
                baer_paths.append(baer_file_path)
                
output_baer_mosaic_path = os.path.join(
    combined_training_data_output_directory, 
    str(path) + "_" + str(row), str(year) + "_baer_mosaic.img")

proc = subprocess.Popen(['C:\Program Files\GDAL\gdal_merge.py', 
                      '-o', output_baer_mosaic_path, 
                      '-pct',
                      '-ul_lr'] + [str(x) for x in extent] + 
                      baer_paths,
                      stdout=subprocess.PIPE, shell=True)
out, err = proc.communicate()

# Mask BAER data with our Data Mask
output_baer_mosaic_masked_path = os.path.join(
     combined_training_data_output_directory, 
     str(path) + "_" + str(row), str(year) + "_baer_mosaic_masked.img")

merged_baer_ds = gdal.Open(output_baer_mosaic_path,0)
merged_baer_band = merged_baer_ds.GetRasterBand(1)
merged_baer_data = merged_baer_band.ReadAsArray(0,0)
    
baer_mask = merged_baer_data * data_mask_data
baer_mask[baer_mask == 1] = 0
baer_mask[baer_mask == 2] = 0
baer_mask[data_mask_data == 0] = 255

output_ds = driver.CreateCopy(output_baer_mosaic_masked_path, 
                              merged_baer_ds, 0)
output_ds.GetRasterBand(1).WriteArray(baer_mask)
output_ds.GetRasterBand(1).SetNoDataValue(255)


baer_mask

detection_data_ds = gdal.Open(destination_path,0)
detection_data = merged_baer_ds.GetRasterBand(1).ReadAsArray(0,0)

print baer_mask[baer_mask == 1] * baer_mask[detection_data == 1]











output_ds = None
    
    
        
#proc = subprocess.Popen(['C:\Program Files\GDAL\gdal_calc.py', 
#                         '-A', output_baer_mosaic_path, 
#                         '-B', data_mask_path, 
#                         '--calc', 'A*B', 
#                             '--outfile=' + output_baer_mosaic_masked_path, 
#                             '--overwrite'], 
#                             stdout=subprocess.PIPE, shell=True)
#                             
#out, err = proc.communicate()
#    
#detected_disturbance_band = detected_disturbance_ds.GetRasterBand(1)
    

    
    
    
    
    
    
