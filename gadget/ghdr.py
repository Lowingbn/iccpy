from sys import argv
from iccpy.gadget import load_snapshot
from iccpy.gadget.binary_snapshot_io import header_names


def dump(name):
    print '\nReading', name, '\n'
    snap = load_snapshot(name)
    print snap
    print '\n'
    header = snap.header
    
    for hname in header_names[:-1]:
        if hname in ['buffer', 'num_particles']:
            # these are stripped of the object
            continue
        dta = getattr(snap.header, hname)
        if len(dta)==1:
            print '%25s'%(hname,), dta[0]
        else:
            print '%25s'%(hname,), dta


if __name__=='__main__':

    if len(argv)!=2:
        print 'Usage: python ghdr.py SNAPSHOT\n where SNAPSHOT is e.g. snapdir_002/snap_002.0 or ICs/IC128.dat'
        exit(0)
    dump(argv[1])
    
