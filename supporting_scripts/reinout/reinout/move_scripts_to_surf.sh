#!/bin/bash

echo "---START----" 

#sshpass -p 'password' 

scp -r ~/scripts/reinout/*.py lofarvwf-jdejong@spider.surfsara.nl:/home/lofarvwf-jdejong/scripts/

echo "---END----"