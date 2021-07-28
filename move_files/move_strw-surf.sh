#!/bin/bash

#THIS SCRIPT IS WRITTEN FOR STRW SERVERS
#With this script you can transfer the files to extract from strw to surf
echo "---START----" 

FROM=$1 #/disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
TO=$2 #/project/lofarvwf/Share/jdejong/data/L626678

cd /net/tussenrijn/data2/jurjendejong/

tar -C ${FROM} -cvzf data_archive.tar.gz --absolute-names image_full_ampphase_di_m.NS.tessel.reg image_full_ampphase_di_m.NS.mask01.fits image_full_ampphase_di_m.NS.DicoModel image_dirin_SSD_m.npy.ClusterCat.npy DDS3_full_*.01_merged.npz DDS3_full_*.01_smoothed.npz *_uv.pre-cal_*.pre-cal.ms.archive SOLSDIR

scp -r /net/tussenrijn/data2/jurjendejong/L626678/data_archive.tar.gz lofarvwf-jdejong@spider.surfsara.nl:${TO}

echo "---END----" 