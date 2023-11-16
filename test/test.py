import os
from glob import glob

# os.system('source /home/jurjen/Documents/Python/lofar_helpers/test/test.sh')

os.system('rm *.txt')

while len(glob('.txt'))<5:
    os.system('source /home/jurjen/Documents/Python/lofar_helpers/test/test.sh')