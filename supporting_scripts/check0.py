import tables

def has0coordinates(h5):
    try:
        h5 = tables.open_file(h5)
        for c in h5.root.sol000.source[:]:
            x, y = c[0][1][0], c[0][1][1]
            if x==0. and y==0.:
                h5.close()
                return True
        h5.close()
    except:
        pass
    return False