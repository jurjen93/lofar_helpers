#!/bin/bash

#THIS SCRIPT IS WRITTEN FOR STRW SERVERS
#FOR WHEN EXTRACT IS RUN ON STRW

FROM=$1 #/net/tussenrijn/data2/jurjendejong/L626678/extract
TO=$2 #/project/lofarvwf/Share/jdejong/output/L626678/extract

scp -r ${FROM}/boxes ${TO}
scp -r ${FROM}/*box_*.dysco.sub.shift.avg.weights.ms.archive0 lofarvwf-jdejong@spider.surfsara.nl:${TO}