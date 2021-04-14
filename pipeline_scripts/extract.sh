#!/bin/bash
#SBATCH -c 32

echo "----------START----------" 

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts

SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong

#start box number
START_N=1

#REORGANIZE FILES
cp ~/scripts/lofar_helpers/h5_merger.py ~/scripts

#MAKE DIRECTORIES
if test -f "${TO}"; then
  echo "${TO} exists"
else
  mkdir ${TO}
  echo "Created ${TO}"
mkdir ${TO}/extract
echo "Created ${TO}/extract"

#MOVE NEEDED FILES
echo "Moving files to ${FOLDER}/extract and untar..."
cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/data_archive.tar.gz ${TO}/extract
rm -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/data_archive.tar.gz

echo "Succesfully finished moving files..."

cd ${TO}/extract
tar -zxvf data_archive.tar.gz
echo"Untarred succesfuly..."

#CREATE BOXES
echo "Create boxes..."
python3 ${SCRIPT_FOLDER}/lofar_helpers/make_boxes.py -f ${TO}/extract/image_full_ampphase_di_m.NS.app.restored.fits -l ${TO} -i false
echo "Succesfully created boxes..."

END_N=$(ls -dq ${TO}/boxes/box*.reg | wc -l)
END_N=1

#EXTRACT
echo "-----STARTED EXTRACT-----"
for ((i=${START_N};i<=${END_N};i++)); do 
sinularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/sub-sources-outside-region.py -b ${TO}/boxes/box_${i}.reg --overwriteoutput -p box_${i} 
echo "Extracted box_${i}"
done
echo "-----FINISHED EXTRACT-----"

echo "----------END----------"