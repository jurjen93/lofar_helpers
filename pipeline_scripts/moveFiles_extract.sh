#!/bin/bash

echo "---START----" 

FROM=$1 #/disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
TO=$2 #lofarvwf-jdejong@spider.surfsara.nl:/project/lofarvwf/Share/jdejong/data/L626678

#sshpass -p 'password' 
# scp /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/sub-sources-outside-region.py lofarvwf-jdejong@spider.surfsara.nl:/home/lofarvwf-jdejong/scripts
# scp /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/runwscleanLBauto.py lofarvwf-jdejong@spider.surfsara.nl:/home/lofarvwf-jdejong/scripts
# scp /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/runwscleanLBautoR.py lofarvwf-jdejong@spider.surfsara.nl:/home/lofarvwf-jdejong/scripts

cd /net/tussenrijn/data2/jurjendejong/L626678

tar -cvzf data_archive.tar.gz --absolute-names ${FROM}/image_full_ampphase_di_m.NS.tessel.reg ${FROM}/image_full_ampphase_di_m.NS.mask01.fits ${FROM}/image_full_ampphase_di_m.NS.DicoModel ${FROM}/image_dirin_SSD_m.npy.ClusterCat.npy ${FROM}/DDS3_full_*.01_merged.npz ${FROM}/DDS3_full_*.01_smoothed.npz ${FROM}/*_uv.pre-cal_*.pre-cal.ms.archive ${FROM}/SOLSDIR ${TO}

scp -r /net/tussenrijn/data2/jurjendejong/L626678/data_archive.tar.gz ${TO}

echo "---END----" 