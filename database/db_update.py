"""
Check database if box needs to be done or not

Script example to put selfcal in progress status for box 21:
python database.py --box=21 --selfcal_inprogress

Use following status:
* DONE --> for fully finished
* INPROGRESS --> if selfcal or extract is in progress
* None --> if not existing

CHECK https://github.com/mhardcastle/lotss-query/blob/master/surveys_db.py FOR FULL CODE
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

from recalibration_db import SurveysDB
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--box', type=str, help='box number', required=True)
parser.add_argument('--check_selfcal', action='store_true', help='check database status for selfcal status')
parser.add_argument('--check_extract', action='store_true', help='check database status for extracted status')
parser.add_argument('--selfcal_inprogress', action='store_true', help='update selfcal status to INPROGRESS')
parser.add_argument('--extract_inprogress', action='store_true', help='update extract status to INPROGRESS')
parser.add_argument('--selfcal_done', action='store_true', help='update selfcal status to DONE')
parser.add_argument('--extract_done', action='store_true', help='update extract status to DONE')
args = parser.parse_args()

TABLE = 'recalibrating'
BOX = 'box_'+str(args.box)

sdb = SurveysDB()

# check if box exists in database
r = sdb.db_get(TABLE, BOX)
if r:
    print(BOX+' exists in recalibration database')
    # check if extracted and self-calibrated
    if args.check_extract:
        if r['extract_status']:
            print(BOX + ': extract status --> '+r['extract_status'])
        else:
            print(BOX + ': is not extracted')
    if args.check_selfcal:
        if r['selfcal_status']:
            print(BOX + ': selfcal status --> '+r['selfcal_status'])
        else:
            print(BOX + ': is not self-calibrated')
else:
    query = 'INSERT INTO recalibrating(id) values '+BOX
    print('EXECUTING QUERY:\n' + query)
    sdb.execute(query)

# update status
sdb.readonly = False
if args.selfcal_inprogress or args.extract_inprogress or args.selfcal_done or args.extract_done:
    if args.selfcal_inprogress:
        query = 'UPDATE recalibrating SET selfcal_status=INPROGRESS WHERE id='+BOX
    if args.extract_inprogress:
        query = 'UPDATE recalibrating SET extract_status=INPROGRESS WHERE id='+BOX
    if args.selfcal_done:
        query='UPDATE recalibrating SET selfcal_status=DONE WHERE id='+BOX
    if args.extract_done:
        query='UPDATE recalibrating SET extract_status=DONE WHERE id='+BOX
    print('EXECUTING QUERY:\n'+query)
    sdb.execute(query)

sdb.close()