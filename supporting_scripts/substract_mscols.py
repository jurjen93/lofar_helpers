import pyrap.tables as pt
import argparse

parser = argparse.ArgumentParser(
    description='Make images from extraction run. Requires working version of the DR2-pipeline software and WSClean (Oct 2018 or newer')
parser.add_argument('--ms', nargs='*', help='msfile(s)')
parser.add_argument('--colname', help='col name in ms')
args = parser.parse_args()

outcolumn = args.colname

for ms in args.ms:
    ts = pt.table(ms, readonly=False)
    colnames = ts.colnames()
    if outcolumn not in colnames:
        desc = ts.getcoldesc('DATA')
        desc['name'] = outcolumn
        ts.addcols(desc)
        ts.close()  # to write results

    else:
        print(outcolumn, ' already exists')
        ts.close()

for ms in args.ms:
    ts = pt.table(ms, readonly=False)
    colnames = ts.colnames()
    if 'CORRECTED_DATA' in colnames:
        data = ts.getcol('CORRECTED_DATA')
    else:
        data = ts.getcol('DATA')
    model = ts.getcol('MODEL_DATA')
    ts.putcol(outcolumn, data - model)
    ts.close()
