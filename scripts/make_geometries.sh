scriptsdir=`pwd`
maindir=${scriptsdir#scripts}

echo $maindir
cd $maindir

mkdir -p geometries
cp templates/*.inp geometries

cd geometries

generate_grid_points.py > generate_grid_settings.out
