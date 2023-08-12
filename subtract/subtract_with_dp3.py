import numpy as np
import sys
import os
import casacore.tables as ct
import tables
import re
import pandas as pd
from subprocess import check_output, STDOUT

def get_largest_divider(inp, max=1000):
    """
    Get largest divider

    :param inp: input number
    :param max: max divider

    :return: largest divider from inp bound by max
    """
    for r in range(max)[::-1]:
        if inp % r == 0:
            return r
    sys.exit("ERROR: code should not arrive here.")

def isfloat(num):
    """
    Check if value is a float
    """
    try:
        float(num)
        return True
    except ValueError:
        return False

def parse_history(ms, hist_item):
    """
    Grep specific history item from MS

    :param ms: measurement set
    :param hist_item: history item

    :return: parsed string
    """
    hist = os.popen('taql "SELECT * FROM '+ms+'::HISTORY" | grep '+hist_item).read().split(' ')
    for item in hist:
        if hist_item in item and len(hist_item)<=len(item):
            return item
    print('WARNING:' + hist_item + ' not found')
    return None


def get_time_preavg_factor(ms: str = None):
    """
    Get time pre-averaging factor (given by demixer.timestep)

    :param ms: measurement set

    :return: averaging integer
    """
    parse_str = "demixer.timestep="
    parsed_history = parse_history(ms, parse_str)
    avg_num = re.findall(r'\d+', parsed_history.replace(parse_str, ''))[0]
    if avg_num.isdigit():
        factor = int(float(avg_num))
        if factor!=1:
            print("WARNING: " + ms + " time has been pre-averaged with factor "+str(factor)+". This might cause time smearing effects.")
        return factor
    elif isfloat(avg_num):
        factor = float(avg_num)
        print("WARNING: parsed factor in " + ms + " is not a digit but a float")
        return factor
    else:
        print("WARNING: parsed factor in " + ms + " is not a float or digit")
        return None

class SubtractDP3:

    def __init__(self, mslist: list = None):
        self.mslist = mslist
        self.cmd = ['DP3',
                   'msin.missingdata=True',
                   'msin.orderms=False',
                   'msout.storagemanager=dysco']
        self.steps = []

    @staticmethod
    def isfulljones(h5: str = None):
        """
        Verify if file is fulljones

        :param h5: h5 file
        """
        T = tables.open_file(h5)
        soltab = list(T.root.sol000._v_groups.keys())[0]
        if 'pol' in T.root.sol000._f_get_child(soltab).val.attrs["AXES"].decode('utf8'):
            if T.root.sol000._f_get_child(soltab).pol[:].shape[0] == 4:
                T.close()
                return True
        T.close()
        return False

    def make_template_modelcolumn(self):
        """
        Make template model column with 0 values
        """

        for ms in self.mslist:

            ts = ct.table(ms, readonly=False)
            colnames = ts.colnames()

            if "MODEL_DATA" not in colnames:
                # get column description from DATA
                desc = ts.getcoldesc('DATA')
                # create output column
                desc['name'] = "MODEL_DATA"
                # create template for output column
                ts.addcols(desc)

            else:
                print("WARNING: MODEL_DATA already exists")
                # get number of rows
                nrows = ts.nrows()
                # make sure every slice has the same size
                best_slice = get_largest_divider(nrows, 1000)
                for c in range(0, nrows, best_slice):
                    model = ts.getcol('MODEL_DATA', startrow=c, nrow=best_slice)
                    ts.putcol('MODEL_DATA', model*0, startrow=c, nrow=best_slice)

    def predict(self,
                sourcedb: list = None,
                h5parm: list = None):
        """
        Predict with DP3 (see https://dp3.readthedocs.io/en/latest/steps/Predict.html)

        :param sourcedb: sky model
        :param h5parm: h5 solutions
        """


        for n, source in enumerate(sourcedb):

            self.steps.append(f'beam{n}')
            self.steps.append(f'predict{n}')

            pnum = source.split('/')[-1].split("_")[0]
            print(pnum, h5parm)
            h5 = [h5 for h5 in h5parm if pnum in h5][0]

            H = tables.open_file(h5)
            direction = H.root.sol000.source[:]['dir'] % (2*np.pi)
            direction *= 360/(2*np.pi)

            self.cmd+= [f'predict{n}.type=predict',
                        f'predict{n}.sourcedb={source}',
                        f'predict{n}.applycal.steps=[amp, phase]',
                        f'predict{n}.applycal.amp.correction=amplitude000',
                        f'predict{n}.applycal.phase.correction=phase000',
                        f'predict{n}.applycal.parmdb={h5}',
                        f'predict{n}.operation=add',
                        f'beam{n}.type=applybeam',
                        f'beam{n}.directions=[{direction[0]}deg,{direction[1]}deg]'
                        'msout.datacolumn=MODEL_DATA',
                        'msin.datacolumn=MODEL_DATA',
                        'msout=.']

        self.cmd += ['steps=' + str(self.steps).replace(" ", "").replace("\'", "")]
        print('\n'.join(self.cmd))

        return self

    def subtract_col(self, out_column: str = None):

        """
        Subtract column in Measurement Set
        :param out_column: out column name
        """

        for ms in self.mslist:
            print('Subtract ' + ms)
            ts = ct.table(ms, readonly=False)
            colnames = ts.colnames()

            if "MODEL_DATA" not in colnames:
                sys.exit(f"ERROR: MODEL_DATA does not exist in {ms}.\nThis is most likely due to a failed predict step.")

            if out_column not in colnames:
                # get column description from DATA
                desc = ts.getcoldesc('DATA')
                # create output column
                desc['name'] = out_column
                # create template for output column
                ts.addcols(desc)

            else:
                print(out_column, ' already exists')

            # get number of rows
            nrows = ts.nrows()
            # make sure every slice has the same size
            best_slice = get_largest_divider(nrows, 1000)
            for c in range(0, nrows, best_slice):
                if c==0:
                    print('SUBTRACT --> DATA - MODEL_DATA')
                data = ts.getcol('DATA', startrow=c, nrow=best_slice)
                model = ts.getcol('MODEL_DATA', startrow=c, nrow=best_slice)
                ts.putcol(out_column, data - model, startrow=c, nrow=best_slice)
            ts.close()

        return self

    def moreDP3(self,
                phaseshift: str = None,
                freqavg: str = None,
                timeavg: str = None, concat: bool = None,
                applybeam: bool = None,
                applycal_h5: str = None):

        """
        Run DP3 command

        :param phaseshift: do phase shift to specific center
        :param freqavg: frequency averaging
        :param timeavg: time averaging
        :param concat: concat the measurement sets
        :param applybeam: apply beam in phaseshifted phase center (or otherwise center of field)
        :param applycal_h5: applycal solution file
        """

        self.cmd += ['msin.datacolumn=SUBTRACT_DATA']

        # 1) PHASESHIFT
        if phaseshift is not None:
            phasecenter = phaseshift.replace('[', '').replace(']', '').split(',')
            phasecenter = f'[{phasecenter[0]},{phasecenter[1]}]'
            self.steps.append('ps')
            self.cmd += ['ps.type=phaseshifter',
                        'ps.phasecenter=' + phasecenter]

        # 2) APPLY BEAM
        if applybeam:
            self.steps.append('beam')
            self.cmd += ['beam.type=applybeam',
                        'beam.direction=[]',
                        'beam.updateweights=True']

        # 3) APPLYCAL
        if applycal_h5 is not None:
            # add fulljones solutions apply
            if self.isfulljones(applycal_h5):
                self.steps.append('ac')
                self.cmd += ['ac.type=applycal',
                            'ac.parmdb=' + applycal_h5,
                            'ac.correction=fulljones',
                            'ac.soltab=[amplitude000,phase000]']
            # add non-fulljones solutions apply
            else:
                ac_count = 0
                T = tables.open_file(applycal_h5)
                for corr in T.root.sol000._v_groups.keys():
                    self.cmd += [f'ac{ac_count}.type=applycal',
                                f'ac{ac_count}.parmdb={applycal_h5}',
                                f'ac{ac_count}.correction={corr}']
                    self.steps.append(f'ac{ac_count}')
                    ac_count += 1

        # 4) AVERAGING
        if freqavg is not None or timeavg is not None:
            self.steps.append('avg')
            self.cmd += ['avg.type=averager']
            if freqavg is not None:
                if str(freqavg).isdigit() or not str(freqavg)[-1].isalpha():
                    self.cmd += [f'avg.freqstep={int(freqavg)}']
                else:
                    self.cmd += [f'avg.freqresolution={freqavg}']
            if timeavg is not None:
                if str(timeavg).isdigit():
                    self.cmd += [f'avg.timestep={int(timeavg)}']
                else:
                    self.cmd += [f'avg.timeresolution={timeavg}']

        self.cmd += ['steps=' + str(self.steps).replace(" ", "").replace("\'", "")]

        if concat:
            self.cmd += [f'msin={",".join(self.mslist)}',
                        'msout=subtract_concat.ms']


        else:
            for n, ms in enumerate(self.mslist):
                self.cmd+=[f'msin={ms}', f'msout=sub_{ms}']

        print('\n'.join(self.cmd))

        return self

    def run(self, type=''):
        """
        Run DP3 command

        :param type: type name
        """

        for n, ms in enumerate(self.mslist):
            self.cmd += [f'msin={ms}', f'msout=sub_{ms}']
            dp3_cmd = open(f"dp3{type}_{n}.cmd", "w")
            dp3_cmd.write('\n'.join(self.cmd))
            dp3_cmd.close()
            check_output(' '.join(self.cmd), stderr=STDOUT, shell=True)

        self.cmd = ['DP3',
                   'msin.missingdata=True',
                   'msin.orderms=False',
                   'msout.storagemanager=dysco']

        return self



if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Subtract region with WSClean')
    parser.add_argument('--mslist', nargs='+', help='measurement sets', required=True)
    parser.add_argument('--sourcedb', nargs='+', help='source models', required=True)
    parser.add_argument('--region', type=str, help='region file', required=True)
    parser.add_argument('--output_name', type=str, help='name of output files (default is model image name)')
    parser.add_argument('--skip_predict', action='store_true',
                        help='skip predict and do only subtract')
    parser.add_argument('--h5parm_predict', nargs='+', help='h5 solution files corresponding with sourcedb')
    parser.add_argument('--phasecenter', type=str,
                        help='phaseshift to given point (example: --phaseshift 16h06m07.61855,55d21m35.4166)')
    parser.add_argument('--freqavg', type=str, help='frequency averaging')
    parser.add_argument('--timeavg', type=str, help='time averaging')
    parser.add_argument('--concat', action='store_true', help='concat MS')
    parser.add_argument('--applybeam', action='store_true', help='apply beam in phaseshift center or center of field')
    parser.add_argument('--applycal', action='store_true', help='applycal after subtraction and phaseshifting')
    parser.add_argument('--applycal_h5', type=str, help='applycal solution file')
    parser.add_argument('--print_only_commands', action='store_true', help='only print commands for testing purposes')
    parser.add_argument('--forwidefield', action='store_true',
                        help='will search for the polygon_info.csv file to extract information from')
    args = parser.parse_args()

    Subtract = SubtractDP3(args.mslist)

    # --forwidefield --> will read averaging and phasecenter from polygon_info.csv
    if args.forwidefield:
        if os.path.isfile('polygon_info.csv'):
            polygon_info = pd.read_csv('polygon_info.csv')
        elif os.path.isfile('../polygon_info.csv'):
            polygon_info = pd.read_csv('../polygon_info.csv')
        elif os.path.isfile('../../polygon_info.csv'):
            polygon_info = pd.read_csv('../../polygon_info.csv')
        elif os.path.isfile('../../../polygon_info.csv'):
            polygon_info = pd.read_csv('../../../polygon_info.csv')
        else:
            sys.exit('ERROR: using --forwidefield option needs polygon_info.csv file to read polygon information from')

        t = ct.table(args.mslist[0] + "::SPECTRAL_WINDOW")
        channum = len(t.getcol("CHAN_FREQ")[0])
        t.close()

        polygon = polygon_info.loc[polygon_info.polygon_file == args.region.split('/')[-1]]
        try:
            phasecenter = polygon['poly_center'].values[0]
        except AttributeError:
            print('WARNING: no poly center in polygon_info.csv, use dir instead.')
            phasecenter = polygon['dir'].values[0]
        except KeyError:
            print('WARNING: no poly center in polygon_info.csv, use dir instead.')
            phasecenter = polygon['dir'].values[0]

        # take only averaging factors that are channum%avg==0
        avg = get_largest_divider(channum, int(polygon['avg'].values[0]))

        freqavg = int(avg)
        try:
            # if there is pre averaging done on the ms, we need to take this into account
            timeavg = int(freqavg/get_time_preavg_factor(args.mslist[0]))
        except:
            timeavg = int(freqavg)
        dirname = polygon['dir_name'].values[0]

    else:
        phasecenter = args.phasecenter
        freqavg = args.freqavg
        timeavg = args.timeavg
        dirname = None

    if not args.skip_predict:
        print('############## PREDICT ##############')
        Subtract.make_template_modelcolumn()
        Subtract.predict(sourcedb=args.sourcedb, h5parm=args.h5parm_predict)
        if not args.print_only_commands:
            Subtract.run(type='predict')
            Subtract.subtract_col('SUBTRACT_DATA')

    if args.phasecenter is not None or \
            args.freqavg is not None or \
            args.timeavg is not None or \
            args.concat is not None or \
            args.applybeam is not None or \
            args.applycal is not None:
        print('############## RUN DP3 ##############')
        if args.applycal_h5 is not None:
            applycalh5 = args.applycal_h5
        elif args.applycal and not args.applycal_h5:
            sys.exit("ERROR: need a solution file for applycal (give with --applycal_h5)")
        else:
            applycalh5 = None

        Subtract.moreDP3(phaseshift=phasecenter, freqavg=freqavg, timeavg=timeavg,
                       concat=args.concat, applybeam=args.applybeam, applycal_h5=applycalh5)
        if not args.print_only_commands:
            Subtract.run(type='phaseshift')
