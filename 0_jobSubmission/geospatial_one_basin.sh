#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%j.out
#==================================================================================
# CAMELS-SPAT geospatial data for ONE basin on ARC: cut the North America datasets
# to the basin outline. Submit from this folder, one job per basin:
#     sbatch -J geo_MX_22008 geospatial_one_basin.sh MX_22008
# Output: $OUT/basin_data/<basin>/geospatial/<dataset>/...   Log: logs/geo_<basin>_<jobid>.out
#==================================================================================
set -euo pipefail
basinId=$1

DATA=/work/comphyd_lab/data/_to-be-moved/NorthAmerica_geospatial   # CAMELS-SPAT datasets (Wouter)
V01=/work/comphyd_lab/share/dataset_develop/camels_mx/basin_delineation_v01
OUT=/work/comphyd_lab/share/dataset_develop/camels_mx/camels_spat_data
ENV=/work/comphyd_lab/users/nvasquez/envs/nwam
export PATH=$ENV/bin:$PATH PYTHONNOUSERSITE=1 PROJ_DATA=$ENV/share/proj
unset PROJ_LIB LD_LIBRARY_PATH
umask 002

# Basin shapefile from v01: the nested export if the basin has one, otherwise the normal export
SHP=$V01/$basinId/cwarhm_exports/tdx_nested_preview/shapefiles/catchment/${basinId}_catchment.shp
[ -f $SHP ] || SHP=$V01/$basinId/cwarhm_exports/tdx_clipped/shapefiles/catchment/${basinId}_catchment.shp
[ -f $SHP ] || { echo "No v01 shapefile for $basinId (export missing or failed)"; exit 1; }

echo "== $basinId  $(date)"
echo "   shapefile: $SHP"
cd $(dirname $SLURM_SUBMIT_DIR)/8_geospatial_data
python 2_subset_geospatial_one_basin.py $basinId $SHP $DATA $OUT
echo "== $basinId DONE $(date): $(du -sh $OUT/basin_data/$basinId | cut -f1), $(find $OUT/basin_data/$basinId -type f | wc -l) files"
