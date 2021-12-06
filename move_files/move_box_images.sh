#!/bin/bash

for ((N=1;N<=71;N++))
do
  scp lofarvwf-jdejong@spider.surfsara.nl:/project/lofarvwf/Share/jdejong/output/A399/selfcal/box_${N}/image_007.app.restored.fits /net/tussenrijn/data2/jurjendejong/A399/box_images/box_${N}.fits
done

