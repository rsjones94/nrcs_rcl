"""
The RCL (Riparian Classification from LidAR) is an Arcpy tool that creates a landcover classification raster using an
ESRI LAS dataset (.lasd) as its input.
"""

import os
from os import listdir
from os.path import isfile, join

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


inlas = arcpy.GetParameterAsText(0) # folder of las files
out_folder = arcpy.GetParameterAsText(1) # folder to create and write to
class_type = arcpy.GetParameterAsText(2) # binary or ternary
z_factor = arcpy.GetParameter(3) # elevation adjustment factor

if class_type not in ['binary', 'ternary']:
    raise Exception('Invalid classification scheme')


######################



progress = 0

support_folder = os.path.join(out_folder,'support')
os.mkdir(out_folder)
os.mkdir(support_folder)

arcpy.AddMessage("Generating footprint and .lasd")
# make footprint
footprint_name = os.path.join(support_folder, 'las_footprint.shp')
files = [f for f in listdir(inlas) if isfile(join(inlas, f))]
spatial_ref = arcpy.Describe(os.path.join(inlas,files[0])).spatialReference
arcpy.PointFileInformation_3d(inlas, footprint_name, 'LAS', '.las', spatial_ref)

# make lasd
inlasd = os.path.join(support_folder, 'support.lasd')
arcpy.CreateLasDataset_management(inlas, inlasd)

ground_files = {'multipoint':os.path.join(support_folder,'ground_multipoint.shp'),
                'tin':os.path.join(support_folder,'ground_tin.adf'),
                'raster':os.path.join(support_folder,'dem.tif'),
                'rawraster':os.path.join(support_folder,'rawdem.tif'),
                'name':'ground'}
surface_files = {'multipoint':os.path.join(support_folder,'surf_multipoint.shp'),
                'tin':os.path.join(support_folder,'surf_tin.adf'),
                'raster':os.path.join(support_folder,'dsm.tif'),
                'rawraster':os.path.join(support_folder,'rawdsm.tif'),
                'name':'surface'}
height_files = {'multipoint':None,
                'tin':None,
                'raster':os.path.join(support_folder,'dhm.tif'),
                'name':'height'}


arcpy.AddMessage("Generating DSM")
surfaceLyr = arcpy.CreateUniqueName('First Return Layer')
arcpy.MakeLasDatasetLayer_management(in_las_dataset=inlasd,
                                     out_layer=surfaceLyr,
                                     class_code=None,
                                     return_values=['Single Return', 'First of Many'])
arcpy.LasDatasetToRaster_conversion(in_las_dataset=surfaceLyr,
                                    out_raster=surface_files['rawraster'],
                                    value_field='ELEVATION',
                                    interpolation_type="TRIANGULATION NATURAL_NEIGHBOR NO_THINNING CLOSEST_TO_MEAN 0", #'TRIANGULATION Linear {point_thinning_type} {point_selection_method} {resolution}',
                                    data_type='FLOAT',
                                    sampling_type='CELLSIZE',
                                    sampling_value=1,
                                    z_factor=z_factor)
arcpy.Clip_management(surface_files['rawraster'], "#", surface_files['raster'],footprint_name, "0", "ClippingGeometry")
arcpy.Delete_management(surface_files['rawraster'])

arcpy.AddMessage("Generating DEM")
groundLyr = arcpy.CreateUniqueName('Last Return Layer')
arcpy.MakeLasDatasetLayer_management(in_las_dataset=inlasd,
                                     out_layer=groundLyr,
                                     class_code=None,
                                     return_values=['Single Return', 'Last of Many'])
arcpy.LasDatasetToRaster_conversion(in_las_dataset=groundLyr,
                                    out_raster=ground_files['rawraster'],
                                    value_field='ELEVATION',
                                    interpolation_type="TRIANGULATION NATURAL_NEIGHBOR NO_THINNING MAXIMUM 0", #'TRIANGULATION Linear {point_thinning_type} {point_selection_method} {resolution}',
                                    data_type='FLOAT',
                                    sampling_type='CELLSIZE',
                                    sampling_value=1,
                                    z_factor=z_factor)
arcpy.Clip_management(ground_files['rawraster'], "#", ground_files['raster'],footprint_name, "0", "ClippingGeometry")
arcpy.Delete_management(ground_files['rawraster'])

#arcpy.management.Delete(surfaceLyr)
#arcpy.management.Delete(groundLyr)

arcpy.AddMessage("Generating DHM")
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

# Con(param, where_clause="Value<=XXX", in_true_raster_or_constant=YYY, in_false_raster_or_constant=ZZZ)

if class_type == 'ternary':
    other_val = 0
    veg_val = 1
    tree_val = 2
    classified = Con(dhm, where_clause="Value<=3.524",
                     in_true_raster_or_constant=Con(dsm_sl, where_clause="Value<=3.205",
                                                    in_true_raster_or_constant=other_val,
                                                    in_false_raster_or_constant=veg_val),
                     in_false_raster_or_constant=tree_val)
elif class_type == 'binary':
    other_val = 0
    tree_val = 1
    classified = Con(dsm_sl, where_clause="Value<=19.241",
                     in_true_raster_or_constant=tree_val,
                     in_false_raster_or_constant=other_val)

classified_path = os.path.join(out_folder, 'classified.tif')
classified.save(classified_path)

arcpy.AddMessage("Classification complete")
arcpy.ResetProgressor()


