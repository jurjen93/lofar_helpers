from surveys_db import SurveysDB

test=True
if test:
    sdb = SurveysDB()
    target='eFEDSJ0845.2+0327'
    currentdict = sdb.get_reprocessing('%s'%target)
    sdb.close()

# def get_reprocessing(self, id):
#     return self.db_get('reprocessing', id)