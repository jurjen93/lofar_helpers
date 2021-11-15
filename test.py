from argparse import ArgumentParser
import sys

parser = ArgumentParser()
parser.add_argument('--test', type=str, default=None)
args = parser.parse_args()

print(args.test)
if args.test:
    if (args.test.startswith("[") and args.test.endswith("]")):
        l = args.test.replace(' ','').replace('[','').replace(']','').split(',')
        for n, v in enumerate(l):
            if not v.isdigit():
                sys.exit('--filter_directions can only have integers in the list.')
            else:
                l[n]=int(v)
    else:
        sys.exit('--test given but no list format. Please pass a list to --test.')
