scriptsdir=`pwd`
maindir=${scriptsdir#scripts}

echo $maindir
cd $maindir

mkdir -p geometries
cp input/*.inp geometries

cd geometries

generate_grid_points.py > generate_grid_settings.out

cd $maindir
./scripts/get_global_coordinates.py
