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

inlas = arcpy.GetParameterAsText(0)  # folder of las files

out_folder = arcpy.GetParameterAsText(1)  # folder to create and write to

training_file = arcpy.GetParameter(2)  # shapefile of training shapes

z_factor = arcpy.GetParameter(3)  # elevation adjustment factor

network_file = arcpy.GetParameter(4)  # stream network
has_clipping = not network_file.isEmpty

outer_buffer_width = arcpy.GetParameter(5)  # riparian buffer width
if outer_buffer_width:
    has_outer = True
else:
    has_outer = False

"""
inner_buffer_width = arcpy.GetParameter(6)  # tolerance buffer width
if inner_buffer_width:
    has_inner = True
else:
    has_inner = False

if has_clipping and not has_inner:
    inner_buffer_width = 0
"""

if has_clipping and not has_outer:
    raise Exception('Riparian buffer width required if stream network is specified')


######################


progress = 0

support_folder = os.path.join(out_folder, 'support')
os.mkdir(out_folder)
os.mkdir(support_folder)

arcpy.AddMessage("Generating footprint and .lasd")

# make lasd
inlasd = os.path.join(support_folder, 'support.lasd')
arcpy.CreateLasDataset_management(inlas, inlasd)

# in the event that a clipping file is specified, make a new folder to dump the extracted las and change the inlas
# and inlasd variables to point to the extracted data
if has_clipping:
    clipping_file = os.path.join(support_folder, 'riparian_buffer.shp')
    arcpy.Buffer_analysis(in_features=network_file,
                          out_feature_class=clipping_file,
                          buffer_distance_or_field=outer_buffer_width,
                          dissolve_option='ALL')  # str(outer_buffer_width) + " Unknown"

    extraction_folder = os.path.join(support_folder, 'extraction')
    extraction_las_folder = os.path.join(extraction_folder, 'las')
    os.mkdir(extraction_folder)
    os.mkdir(extraction_las_folder)
    arcpy.ExtractLas_3d(in_las_dataset=inlasd,
                        target_folder=extraction_las_folder,
                        extent=None,
                        boundary=clipping_file)

    inlas = extraction_las_folder
    inlasd = os.path.join(extraction_folder, 'support_extracted.lasd')
    arcpy.CreateLasDataset_management(inlas, inlasd)

# make footprint. if a clipping file was specified, the footprint will be of the extracted las
footprint_name = os.path.join(support_folder, 'las_footprint.shp')
files = [f for f in listdir(inlas) if isfile(join(inlas, f))]
spatial_ref = arcpy.Describe(os.path.join(inlas, files[0])).spatialReference
arcpy.PointFileInformation_3d(inlas, footprint_name, 'LAS', '.las', spatial_ref)

ground_files = {'multipoint': os.path.join(support_folder, 'ground_multipoint.shp'),
                'tin': os.path.join(support_folder, 'ground_tin.adf'),
                'raster': os.path.join(support_folder, 'dem.tif'),
                'rawraster': os.path.join(support_folder, 'rawdem.tif'),
                'name': 'ground'}
surface_files = {'multipoint': os.path.join(support_folder, 'surf_multipoint.shp'),
                 'tin': os.path.join(support_folder, 'surf_tin.adf'),
                 'raster': os.path.join(support_folder, 'dsm.tif'),
                 'rawraster': os.path.join(support_folder, 'rawdsm.tif'),
                 'name': 'surface'}
height_files = {'multipoint': None,
                'tin': None,
                'raster': os.path.join(support_folder, 'dhm.tif'),
                'name': 'height'}
intensity_files = {'multipoint': None,
                   'tin': None,
                   'raster': os.path.join(support_folder, 'intensity.tif'),
                   'rawraster': os.path.join(support_folder, 'rawintensity.tif'),
                   'name': 'intensity'}

arcpy.AddMessage("Generating DSM")
surfaceLyr = arcpy.CreateUniqueName('First Return Layer')
arcpy.MakeLasDatasetLayer_management(in_las_dataset=inlasd,
                                     out_layer=surfaceLyr,
                                     class_code=None,
                                     return_values=['Single Return', 'First of Many'])
arcpy.LasDatasetToRaster_conversion(in_las_dataset=surfaceLyr,
                                    out_raster=surface_files['rawraster'],
                                    value_field='ELEVATION',
                                    interpolation_type="TRIANGULATION NATURAL_NEIGHBOR NO_THINNING CLOSEST_TO_MEAN 0",
                                    # 'TRIANGULATION Linear {point_thinning_type} {point_selection_method} {resolution}',
                                    data_type='FLOAT',
                                    sampling_type='CELLSIZE',
                                    sampling_value=1,
                                    z_factor=z_factor)

arcpy.AddMessage("Generating DEM")
groundLyr = arcpy.CreateUniqueName('Last Return Layer')
arcpy.MakeLasDatasetLayer_management(in_las_dataset=inlasd,
                                     out_layer=groundLyr,
                                     class_code=None,
                                     return_values=['Single Return', 'Last of Many'])
arcpy.LasDatasetToRaster_conversion(in_las_dataset=groundLyr,
                                    out_raster=ground_files['rawraster'],
                                    value_field='ELEVATION',
                                    interpolation_type="TRIANGULATION NATURAL_NEIGHBOR NO_THINNING MAXIMUM 0",
                                    # 'TRIANGULATION Linear {point_thinning_type} {point_selection_method} {resolution}',
                                    data_type='FLOAT',
                                    sampling_type='CELLSIZE',
                                    sampling_value=1,
                                    z_factor=z_factor)

arcpy.AddMessage("Generating Intensity Raster")
allLyr = arcpy.CreateUniqueName('All Return Layer')
arcpy.MakeLasDatasetLayer_management(in_las_dataset=inlasd,
                                     out_layer=allLyr,
                                     class_code=None,
                                     return_values=None)
arcpy.LasDatasetToRaster_conversion(in_las_dataset=allLyr,
                                    out_raster=intensity_files['rawraster'],
                                    value_field='INTENSITY',
                                    interpolation_type="TRIANGULATION NATURAL_NEIGHBOR NO_THINNING MAXIMUM 0",
                                    # 'TRIANGULATION Linear {point_thinning_type} {point_selection_method} {resolution}',
                                    data_type='FLOAT',
                                    sampling_type='CELLSIZE',
                                    sampling_value=1,
                                    z_factor=z_factor)

"""
if has_clipping:
    merged_intersect = os.path.join(extraction_folder, 'footprint_intersection.shp')
    arcpy.Intersect_analysis(in_features=[footprint_name, clipping_file],
                             out_feature_class=merged_intersect)
    cutter = merged_intersect
else:
    cutter = footprint_name


arcpy.Clip_management(surface_files['rawraster'], "#", surface_files['raster'], cutter, "0", "ClippingGeometry")
arcpy.Clip_management(ground_files['rawraster'], "#", ground_files['raster'], cutter, "0", "ClippingGeometry")
arcpy.Delete_management(surface_files['rawraster'])
arcpy.Delete_management(ground_files['rawraster'])
"""

# arcpy.management.Delete(surfaceLyr)
# arcpy.management.Delete(groundLyr)

arcpy.AddMessage("Generating DHM")
memory_dhm = arcpy.sa.Raster(surface_files['rawraster']) - arcpy.sa.Raster(ground_files['rawraster'])
memory_dhm.save(height_files['rawraster'])


dem = ground_files['rawraster']
dsm = surface_files['rawraster']
dhm = height_files['rawraster']
inten = intensity_files['rawraster']

dem_sl = os.path.join(support_folder, 'demsl.tif')
dsm_sl = os.path.join(support_folder, 'dsmsl.tif')
dhm_sl = os.path.join(support_folder, 'dhmsl.tif')

arcpy.AddMessage("Generating slope rasters")
progress += 1
for ras, ras_sl in zip([dem, dsm, dhm], [dem_sl, dsm_sl, dhm_sl]):
    memory_sl = arcpy.sa.Slope(in_raster=ras,
                               output_measurement='DEGREE',
                               z_factor=None,  # already accounted for at this pint
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

pass

arcpy.AddMessage("Cleaning classifications")

classified_path_raw = os.path.join(support_folder, 'classified_raw.tif')
classified.save(classified_path_raw)

classified_path_maj = os.path.join(support_folder, 'classified_maj.tif')
maj = arcpy.sa.MajorityFilter(classified_path_raw, "EIGHT", "HALF")
maj.save(classified_path_maj)

classified_path_clean = os.path.join(support_folder, 'classified_clean.tif')
cle = arcpy.sa.BoundaryClean(classified_path_maj, "DESCEND", "TWO_WAY")
cle.save(classified_path_clean)

arcpy.AddMessage("Converting to polygons")

classified_path_poly = os.path.join(out_folder, 'classified_poly.shp')
arcpy.RasterToPolygon_conversion(classified_path_clean, classified_path_poly, "SIMPLIFY")

if has_clipping:
    arcpy.AddMessage("Calculating distances")
    arcpy.Near_analysis(classified_path_poly, network_file, outer_buffer_width)

arcpy.AddMessage("Classification complete")
arcpy.ResetProgressor()
