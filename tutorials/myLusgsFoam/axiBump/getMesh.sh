#!/bin/bash
cd ${0%/*} || exit 1                        # Run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions    # Tutorial run functions

# Uncomment the mesh file
mesh=bump_newaxi_91.p3dfmt
#mesh=bump_newaxi_181.p3dfmt
#mesh=bump_newaxi_361.p3dfmt
#mesh=bump_newaxi_721.p3dfmt
#mesh=bump_newaxi_1441.p3dfmt

if [ ! -f ${mesh} ]; then
    wget https://turbmodels.larc.nasa.gov/Axi_Bump_grids/${mesh}.gz
    gzip -d ${mesh}.gz
fi

runApplication plot3dToFoam -noBlank ${mesh}
runApplication transformPoints -yawPitchRoll '(0 0 -90.5)'
runApplication autoPatch -overwrite 60
runApplication createPatch -overwrite
runApplication renumberMesh -overwrite

echo "WARNING: check the boundary condition assignment!"
