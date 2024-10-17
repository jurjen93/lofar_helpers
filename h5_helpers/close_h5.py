import tables
import sys

def list_opened():
    """List current open HDF5 files"""
    h5s = tables.file._open_files.filenames
    return h5s


def force_close_all():
    """Close all currently open HDF5 files by force"""
    if len(tables.file._open_files._handlers) > 0:
        tables.file._open_files.close_all()
    else:
        sys.stderr.write('No open HDF5 files'+'\n')


def force_close(h5):
    """Close indivdual HDF5 file by force"""
    h5s = list(tables.file._open_files._handlers)
    for h in h5s:
        if h.filename == h5:
            sys.stdout.write('Closed --> ' + h5+'\n')
            h.close()
            return
    sys.stderr.write(h5 + ' not found\n')


if __name__ == '__main__':
    force_close_all()