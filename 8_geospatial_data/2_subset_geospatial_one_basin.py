"""Cut the CAMELS-SPAT geospatial datasets to ONE basin (CAMELS-MX).

Same result as 2a (rasters) + 2b (HydroLAKES, GLHYMPS) + 3b (soil porosity and
conductivity), using the same CAMELS-SPAT functions, but it takes the basin
shapefile directly: no metadata CSV and no config file.

    python 2_subset_geospatial_one_basin.py BASIN_ID SHAPEFILE DATA_DIR OUT_DIR

    BASIN_ID   e.g. MX_22008
    SHAPEFILE  basin shapefile (its HRUs are dissolved into one outline)
    DATA_DIR   folder with the North America datasets (CAMELS-SPAT geospatial_temp)
    OUT_DIR    writes OUT_DIR/basin_data/BASIN_ID/shapefiles/lumped/ and .../geospatial/

Files that already exist are skipped, so a failed job can simply be resubmitted.
"""
import os
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio

sys.path.append(str(Path(__file__).resolve().parents[1]))  # repo root, for python_cs_functions
import python_cs_functions as cs

# The 11 datasets CAMELS-SPAT uses (folder names inside DATA_DIR)
DATASETS = ['worldclim', 'lai', 'forest_height', 'glclu2019', 'modis_land', 'lgrip30',
            'merit', 'pelletier', 'soilgrids', 'hydrolakes', 'glhymps']

basin, shapefile, data_dir, out_dir = sys.argv[1:5]
data_dir = Path(data_dir)
basin_dir = Path(out_dir) / 'basin_data' / basin

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

# 2. Cut every dataset to the basin (2a for .tif, 2b for .shp), keeping the folder structure
for ds in DATASETS:
    n = 0
    for root, _, files in os.walk(data_dir / ds):
        for f in sorted(files):
            src = Path(root) / f
            des = basin_dir / 'geospatial' / src.parent.relative_to(data_dir)
            if (des / f).exists():
                continue
            if f.endswith('.tif'):
                cs.subset_geotiff_to_shapefile(src, lump_shp, des, buffer=True)
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
    print(f'  {ds}: {n} files cut', flush=True)

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
