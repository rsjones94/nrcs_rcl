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

from arcpy.sa import Con

'''
inlas = r'F:\gen_model\rcl_testing\las' # folder of las files
out_folder = r'F:\gen_model\rcl_testing\rcl_model' # folder to create and write to
pt_space = 0.6 # point spacing
z_factor = 1 # elevation adjustment factor
#coord_sys = r"Coordinate Systems\Projected Coordinate Systems\UTM\NAD 1983\NAD 1983 UTM Zone 16N.prj"
coord_sys = None
'''

inlas = arcpy.GetParameterAsText(0) # folder of las files
out_folder = arcpy.GetParameterAsText(1) # folder to create and write to
pt_space = arcpy.GetParameter(2) # point spacing
z_factor = arcpy.GetParameter(3) # elevation adjustment factor
coord_sys = arcpy.GetParameter(4)

###

progress = 0

support_folder = os.path.join(out_folder,'support')
os.mkdir(out_folder)
os.mkdir(support_folder)

ground_files = {'multipoint':os.path.join(support_folder,'ground_multipoint.shp'),
                'tin':os.path.join(support_folder,'ground_tin.adf'),
                'raster':os.path.join(support_folder,'dem.tif'),
                'name':'ground'}
surface_files = {'multipoint':os.path.join(support_folder,'surf_multipoint.shp'),
                'tin':os.path.join(support_folder,'surf_tin.adf'),
                'raster':os.path.join(support_folder,'dsm.tif'),
                'name':'surface'}
height_files = {'multipoint':None,
                'tin':None,
                'raster':os.path.join(support_folder,'dhm.tif'),
                'name':'height'}

arcpy.SetProgressorLabel("Generating multipoint files")

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

arcpy.AddMessage("Generating TINs and elevation rasters")
progress += 1

# then we need to create TINs and then make a raster from those
for map in [ground_files, surface_files]:

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
memory_dhm = arcpy.sa.Raster(surface_files['raster']) - arcpy.sa.Raster(ground_files['raster'])
memory_dhm.save(height_files['raster'])


dem = ground_files['raster']
dsm = surface_files['raster']
dhm = height_files['raster']

dem_sl = os.path.join(support_folder,'demsl.tif')
dsm_sl = os.path.join(support_folder,'dsmsl.tif')
dhm_sl = os.path.join(support_folder,'dhmsl.tif')

arcpy.AddMessage("Generating slope rasters")
progress += 1

for ras, ras_sl in zip([dem, dsm, dhm], [dem_sl, dsm_sl, dhm_sl]):
    memory_sl = arcpy.sa.Slope(in_raster=ras,
                               output_measurement='DEGREE',
                               z_factor=None, # already accounted for at this pint
                               method='PLANAR',
                               z_unit=None)
    memory_sl.save(ras_sl)

    # we're replacing the filepaths in the list we're iterating through with in-memory rasters as we go
    # doing this to make the condition statement coming up somewhat easier to parse
    ras = arcpy.sa.Raster(ras)
    ras_sl = arcpy.sa.Raster(ras_sl)

# now that we have all out data we can run the decision tree

arcpy.AddMessage("Classifying cover")
progress += 1

other_val = 3
veg_val = 2
tree_val = 1

classified = Con(dsm_sl, where_clause="Value<=29.117",
                 in_true_raster_or_constant=Con(dsm_sl, where_clause="Value<=2.887",
                                                in_true_raster_or_constant=other_val,
                                                in_false_raster_or_constant=Con(dem_sl, where_clause="Value<=32.663",
                                                                                in_true_raster_or_constant=veg_val,
                                                                                in_false_raster_or_constant=tree_val)),
                 in_false_raster_or_constant=Con(dsm_sl, where_clause="Value<=55.36",
                                                 in_true_raster_or_constant=Con(dem_sl, where_clause="Value<=45.682",
                                                                                in_true_raster_or_constant=veg_val,
                                                                                in_false_raster_or_constant=tree_val),
                                                 in_false_raster_or_constant=tree_val))

classified_path = os.path.join(out_folder, 'classified.tif')
classified.save(classified_path)

arcpy.AddMessage("Classification complete")
arcpy.ResetProgressor()

