#!/bin/bash

#THIS SCRIPT IS WRITTEN FOR STRW SERVERS

FROM=$1 #/net/tussenrijn/data2/jurjendejong/L626678
TO=$2 #lofarvwf-jdejong@spider.surfsara.nl:/project/lofarvwf/Share/jdejong/output/L626678

scp -r ${FROM}/boxes ${TO}
scp -r ${FROM}/extract/*box_*.dysco.sub.shift.avg.weights.ms.archive0 ${TO}/extract/