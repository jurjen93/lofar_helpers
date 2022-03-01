import tables

def list_opened():
    """List current open HDF5 files"""
    h5s = tables.file._open_files.filenames
    print('Currently open --> '+', '.join(h5s))
    return h5s

def force_close_all():
    """Forcefully close all open HDF5 files"""
    if len(tables.file._open_files._handlers)>0:
        tables.file._open_files.close_all()
    else:
        print('No open HDF5 files')

def force_close(h5):
    """Forcefully close indivdual HDF5 file"""
    h5s = list(tables.file._open_files._handlers)
    for h in h5s:
        if h.filename==h5:
            print('Closed by force --> ' + h5)
            h.close()
            return
    print(h5 + ' not found')
