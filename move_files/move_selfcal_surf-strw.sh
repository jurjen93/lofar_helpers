#!/bin/bash

FROM=$1 #/project/lofarvwf/Share/jdejong/output/L626678
TO=$2 #/net/tussenrijn/data2/jurjendejong/L626678

scp lofarvwf-jdejong@spider.surfsara.nl:${FROM}/all_directions.h5 ${TO}