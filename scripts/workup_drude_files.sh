#!/bin/bash
scriptsdir=`pwd`
maindir=${scriptsdir#scripts}
dispersiondir=$maindir/dispersion/
thole_param=2.0

constraints=$maindir/templates/dispersion_base_constraints.index
drude_constraints=$maindir/templates/drude_base_constraints.index
cd $dispersiondir
WD=`pwd`

echo "Getting drude base constraints from the following file:"
echo $drude_constraints
echo "Make sure atomtypes from base constraint files correspond to atomtypes from .clt file!"
echo 

for dir in `ls *clt`
do
    # Get file names and go to dispersion directory.
    NAME=${dir%.clt}
    if [ $NAME == template ]; then continue; fi
    echo $NAME

    path=${PWD#/home/mvanvleet/}
    homedir=$WD
    cd $NAME

    # Check for .pol files, and localize otherwise
    count=`ls -1 *.pol 2>/dev/null | wc -l`
    if [ $count == 0 ]
    then 
        localize.sh $NAME --limit 3 --isotropic --weight 0 --keep
        echo
        echo "(If above line indicates an error, don't worry overmuch)"
    fi 
    cd $homedir


    # Use constraint file to get atomic dispersion coefficients
    #cp $constraints .

    ## mv *.casimir old.casimir 2>/dev/null
    ## main_dispersion $NAME $path/ $(basename $constraints) > ${NAME}_fit_dispersion.out

    ## casimir < ${NAME}*.casimir > ${NAME}.cncoeffs

    # Use drude constraint file to get atomic drude charges
    cp $drude_constraints .
    path=${path#research/working_directory}
    main_drude $NAME $path/ $(basename $drude_constraints) $thole_param > ${NAME}_fit_drude.out

done
