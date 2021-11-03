import tables

def has0coordinates(h5):
    # try:
    h5 = tables.open_file(h5)
    for c in h5.root.sol000.source[:]:
        x, y = c[1]
        if x==0. and y==0.:
            h5.close()
            return True
    h5.close()
    # except:
    #     h5.close()
        # pass
    return False

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--h5', type=str, help='h5 file name')
    args = parser.parse_args()

    if has0coordinates(args.h5):
        print(args.h5+' has 0.0 coordinates')
    else:
        print(args.h5+' has no 0.0 coordinates')