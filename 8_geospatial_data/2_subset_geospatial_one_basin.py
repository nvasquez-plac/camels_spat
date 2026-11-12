"""Cut the CAMELS-SPAT geospatial datasets to ONE basin (CAMELS-MX).

Same result as 2a (rasters) + 2b (HydroLAKES, GLHYMPS) + 3b (soil porosity and
conductivity), using the same CAMELS-SPAT functions, but it takes the basin
shapefile directly: no metadata CSV and no config file.

    python 2_subset_geospatial_one_basin.py BASIN_ID SHAPEFILE DATA_DIR OUT_DIR

    BASIN_ID   e.g. MX_22008
    SHAPEFILE  basin shapefile (its HRUs are dissolved into one outline)
    DATA_DIR   folder with the North America datasets (CAMELS-SPAT geospatial_temp)
    OUT_DIR    writes OUT_DIR/BASIN_ID/shapefiles/lumped/ and OUT_DIR/BASIN_ID/geospatial/

Files that already exist are skipped, so a failed job can simply be resubmitted.
A file that fails is logged and skipped; the job then ends with an error listing them.
"""
import glob
import os
import sys
import time
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from osgeo import gdal
from shapely import make_valid

sys.path.append(str(Path(__file__).resolve().parents[1]))  # repo root, for python_cs_functions
import python_cs_functions as cs

# The 11 datasets CAMELS-SPAT uses (folder names inside DATA_DIR)
DATASETS = ['worldclim', 'lai', 'forest_height', 'glclu2019', 'modis_land', 'lgrip30',
            'merit', 'pelletier', 'soilgrids', 'hydrolakes', 'glhymps']

basin, shapefile, data_dir, out_dir = sys.argv[1:5]
data_dir = Path(data_dir)
basin_dir = Path(out_dir) / basin

# 1. Basin outline = CAMELS-SPAT "lumped" shapefile. EPSG:4326, because the clipping buffers in degrees.
lump_dir = basin_dir / 'shapefiles' / 'lumped'
lump_dir.mkdir(parents=True, exist_ok=True)
lump_shp = lump_dir / f'{basin}_lumped.shp'
hru = gpd.read_file(shapefile)
lump = hru[['geometry']].dissolve().to_crs('EPSG:4326')
lump['basin_id'] = basin
lump['area_km2'] = lump.to_crs('ESRI:102008').area / 1e6  # same equal-area CRS as CAMELS-SPAT
lump.to_file(lump_shp)
print(f'{basin}: {len(hru)} HRUs -> one outline, {lump.area_km2.iloc[0]:.1f} km2', flush=True)

# 2. Cut every dataset to the basin (2a for .tif, 2b for .shp), keeping the folder structure.
#    Rasters: same steps as cs.subset_geotiff_to_shapefile (buffer the outline by half a pixel
#    diagonal, no buffer for MERIT, then gdal.Warp with the outline as cutline), but the buffered
#    outline is made once per pixel size instead of once per file (LAI has 1,984 files).
buffered = {}
def buffered_outline(src):
    if 'merit' in str(src):
        return str(lump_shp)
    ds = gdal.Open(str(src))
    px, py = ds.GetGeoTransform()[1], ds.GetGeoTransform()[5]
    ds = None
    dist = 0.5 * (px**2 + py**2) ** 0.5
    if dist not in buffered:
        tmp = lump_dir / f'{basin}_buffer_{len(buffered)}_TEMP.shp'
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')  # buffering in degrees is intended here
            shp = gpd.read_file(lump_shp)
            shp['geometry'] = shp.buffer(dist).apply(make_valid)
            shp.to_file(tmp)
        buffered[dist] = str(tmp)
    return buffered[dist]

failed = []
for ds_name in DATASETS:
    n, t0 = 0, time.time()
    for root, dirs, files in os.walk(data_dir / ds_name):
        dirs[:] = [d for d in dirs if d != 'raw_reruns']  # leftover SoilGrids download tiles, not used
        for f in sorted(files):
            src = Path(root) / f
            des = basin_dir / 'geospatial' / src.parent.relative_to(data_dir)
            if (des / f).exists():
                continue
            try:
                if f.endswith('.tif'):
                    des.mkdir(parents=True, exist_ok=True)
                    gdal.Warp(destNameOrDestDS=str(des / f), srcDSOrSrcDSTab=str(src),
                              cutlineDSName=buffered_outline(src), cropToCutline=True,
                              copyMetadata=True, dstNodata=None,
                              srcSRS='EPSG:4326', dstSRS='EPSG:4326')
                    n += 1
                elif f.endswith('.shp'):
                    # Read only the basin's bounding box: same result, much faster than all of North America
                    shp = gpd.read_file(src, bbox=tuple(lump.total_bounds))
                    if len(shp) == 0:
                        des.mkdir(parents=True, exist_ok=True)
                        (des / f.replace('.shp', '.txt')).write_text('Source data does not intersect with catchment shapefile')
                    else:
                        cs.subset_shapefile_to_shapefile(shp, src, lump_shp, des)
                    n += 1
            except Exception as e:
                failed.append(f'{src}: {e}')
                print(f'  FAILED {src}: {e}', flush=True)
    print(f'  {ds_name}: {n} files cut in {time.time() - t0:.0f} s', flush=True)

for tmp in buffered.values():  # remove the temporary buffered outlines
    for g in glob.glob(tmp.replace('.shp', '.*')):
        os.remove(g)

# 3. Soil porosity and conductivity from SoilGrids sand and clay (= 3b, Cosby et al. 1984, Table 4)
soil = basin_dir / 'geospatial' / 'soilgrids' / 'raw'
for d in ['0-5cm', '5-15cm', '15-30cm', '30-60cm', '60-100cm', '100-200cm']:
    with rasterio.open(soil / 'sand' / f'sand_{d}_mean.tif') as src:
        sand = np.ma.masked_equal(src.read(), src.nodata)
    with rasterio.open(soil / 'clay' / f'clay_{d}_mean.tif') as src:
        clay = np.ma.masked_equal(src.read(), src.nodata)
        profile = src.profile
        profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=src.nodata)
    sand = sand / 10  # g/kg -> %
    clay = clay / 10
    porosity = (50.5 - 0.142 * sand - 0.037 * clay) / 100                       # fraction
    conductivity = 2.54 * np.power(10, -0.6 + 0.0126 * sand - 0.0064 * clay)  # cm/hr
    for name, values in [('porosity', porosity), ('conductivity', conductivity)]:
        (soil / name).mkdir(parents=True, exist_ok=True)
        with rasterio.open(soil / name / f'{name}_{d}_mean.tif', 'w', **profile) as dst:
            dst.write(values.filled().astype(rasterio.float32))
print('  soil porosity + conductivity: 12 files written', flush=True)

if failed:
    (basin_dir / 'geospatial_failed_files.txt').write_text('\n'.join(failed) + '\n')
    sys.exit(f'{len(failed)} files failed, listed in {basin_dir}/geospatial_failed_files.txt')
