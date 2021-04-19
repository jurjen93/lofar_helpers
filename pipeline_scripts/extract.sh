#!/bin/bash
#SBATCH -c 32

#THIS SCRIPT IS WRITTEN FOR SLURM ON SURFSARA
echo "----------START----------"

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts

SING_IMAGE_1=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_IMAGE_2=/project/lofarvwf/Software/lofar_sksp_fedora27_ddf_slurmfix.sif
SING_BIND=/project/lofarvwf/Share/jdejong

#start box number
START_N=1

#REORGANIZE FILES
cp ~/scripts/lofar_helpers/h5_merger.py ~/scripts
cp ~/scripts/lofar_helpers/merge_selfcals.py ~/scripts

#MOVE NEEDED FILES
#echo "Moving files to ${TO}/extract and untar..."
#cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/data_archive.tar.gz ${TO}/extract
#rm -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/data_archive.tar.gz
#echo "Succesfully finished moving files..."

#cd ${TO}/extract
#tar -zxvf data_archive.tar.gz
#echo"Untarred succesfuly..."

#CREATE BOXES
echo "Create boxes..."
singularity exec -B ${SING_BIND} ${SING_IMAGE_1} python3 ${SCRIPT_FOLDER}/lofar_helpers/make_boxes.py -f ${TO}/extract/image_full_ampphase_di_m.NS.app.restored.fits -l ${TO}
echo "Succesfully created boxes..."

END_N=$(ls -dq ${TO}/boxes/box*.reg | wc -l)
echo "There are ${END_N} boxes"
END_N=1

#EXTRACT
echo "-----STARTED EXTRACT-----"
cd ${TO}/extract
for ((i=${START_N};i<=${END_N};i++)); do 
singularity exec -B ${SING_BIND} ${SING_IMAGE_2} python ${SCRIPT_FOLDER}/sub-sources-outside-region.py -b ${TO}/boxes/box_${i}.reg --overwriteoutput -p box_${i}
echo "Extracted box_${i}"
done
echo "-----FINISHED EXTRACT-----"

echo "----------END----------"