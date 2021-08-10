"""
Check database if box needs to be done or not

Script example to put selfcal in progress status for box 21:
python database.py --box=21 --selfcal_inprogress

Use following status:
* DONE --> for fully finished
* INPROGRESS --> if selfcal or extract is in progress
* None --> if not existing

"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

from surveys_db import SurveysDB
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

BOX = 'box_'+str(args.box)

sdb = SurveysDB()

# check if box exists in database
r = sdb.db_get('recalibrating', BOX)
print(r['selfcal_status'])
print(r['extract_status'])
if r:
    print(BOX+' exists in recalibration database')
    # check if extracted and self-calibrated
    if args.check_extract:
        if r['extract_status']:
            print(BOX+': extract status --> '+r['extract_status'])
        else:
            print(BOX + ': is not extracted')
    if args.check_selfcal:
        if r['selfcal_status']:
            print(BOX+': selfcal status --> '+r['selfcal_status'])
        else:
            print(BOX + ': is not self-calibrated')
else:
    r = sdb.db_create('recalibrating', BOX)
    print('Created new entry: '+BOX)

# update status
r.readonly = False
if args.selfcal_inprogress:
    r['selfcal_status'] = 'INPROGRESS'
    print(BOX + ': selfcal status updated --> ' + r['selfcal_status'])
if args.extract_inprogress:
    r['extract_status'] = 'INPROGRESS'
    print(BOX + ': extract status updated --> ' + r['extract_status'])
if args.selfcal_done:
    r['selfcal_status'] = 'DONE'
    print(BOX + ': selfcal status updated --> ' + r['selfcal_status'])
if args.extract_done:
    r['extract_status'] = 'DONE'
    print(BOX + ': extract status updated --> ' + r['extract_status'])

sdb.close()