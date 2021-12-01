import tables

def has0coordinates(h5):
    h5 = tables.open_file(h5)
    for c in h5.root.sol000.source[:]:
        x, y = c[1]
        if x==0. and y==0.:
            h5.close()
            return True
    h5.close()
    return False

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--h5', type=str, help='h5 file name', default=None)
    parser.add_argument('--folder', type=str, help='folder name with boxes', default=None)
    args = parser.parse_args()

    if args.h5:
        if has0coordinates(args.h5):
            print(args.h5+' has 0.0 coordinates')
        else:
            print(args.h5+' has no 0.0 coordinates')
    if args.folder:
        from glob import glob
        folder = glob(str(args.folder)+'/box_*')
        for box in folder:
            h5files = glob(box+'/*.h5')
            for h5 in h5files:
                if has0coordinates(h5):
                    print(args.h5 + ' has 0.0 coordinates')