"""
The RCL (Riparian Classification from LidAR) is an Arcpy tool that creates a landcover classification raster using an
ESRI LAS dataset (.lasd) as its input.
"""

import os

import arcpy

if arcpy.CheckExtension("3D") == "Available":
    arcpy.CheckOutExtension("3D")
else:
    # Raise a custom exception
    #
    raise Exception('No 3D analyst license available')

if arcpy.CheckExtension("Spatial") == "Available":
    arcpy.CheckOutExtension("Spatial")
else:
    # Raise a custom exception
    #
    raise Exception('No Spatial Analyst license available')


inlas = r'F:\gen_model\rcl_testing\las' # folder of las files
out_folder = r'F:\gen_model\rcl_testing\rcl_model' # folder to create and write to
pt_space = 0.6 # point spacing
z_factor = 1 # elevation adjustment factor
#coord_sys = r"Coordinate Systems\Projected Coordinate Systems\UTM\NAD 1983\NAD 1983 UTM Zone 16N.prj"
coord_sys = None


###

support_folder = os.path.join(out_folder,'support')
os.mkdir(out_folder)
os.mkdir(support_folder)

ground_files = {'multipoint':os.path.join(support_folder,'ground_multipoint.shp'),
                'tin':os.path.join(support_folder,'ground_tin.adf'),
                'raster':os.path.join(support_folder,'dem.tif')}
surface_files = {'multipoint':os.path.join(support_folder,'surf_multipoint.shp'),
                'tin':os.path.join(support_folder,'surf_tin.adf'),
                'raster':os.path.join(support_folder,'dsm.tif')}
height_files = {'multipoint':None,
                'tin':None,
                'raster':os.path.join(support_folder,'dhm.tif')}
name_maps = [ground_files, surface_files]


# we first need to create multipoint files representing the ground (last returns) and surface (first returns)
arcpy.LASToMultipoint_3d(inlas, #input folder
                         ground_files['multipoint'], #output name
                         pt_space, # point spacing
                         None, # class code
                         'LAST_RETURNS', # return number. have to do non keyword parameters because ESRI hilariously made return a parameter keyword which breaks scripts
                         None, # attribute
                         coord_sys, # coordinate system
                         'las', # file suffix
                         z_factor, # z factor
                         False # folder recursion
                         )

arcpy.LASToMultipoint_3d(inlas, #input folder
                         surface_files['multipoint'], #output name
                         pt_space, # point spacing
                         None, # class code
                         '1', # return number. have to do non keyword parameters because ESRI hilariously made return a parameter keyword which breaks scripts
                         None, # attribute
                         coord_sys, # coordinate system
                         'las', # file suffix
                         z_factor, # z factor
                         False # folder recursion
                         )

# then we need to create TINs and then make a raster from those
for map in name_maps:

    in_feats = str(map['multipoint']) + ' Shape.Z masspoints'
    arcpy.CreateTin_3d(out_tin=map['tin'],
                       spatial_reference=coord_sys,
                       in_features=in_feats)
    arcpy.TinRaster_3d(in_tin=map['tin'],
                       out_raster=map['raster'],
                       data_type='FLOAT',
                       method='LINEAR',
                       sample_distance='CELLSIZE 1',
                       z_factor=None # the z factor has already been applied at this point
                       )

# make the height model
dhm = arcpy.sa.Raster(surface_files['raster']) - arcpy.sa.Raster(ground_files['raster'])
dhm.save(height_files['raster'])
