import h5py

class Struct:
    def __init__(self, *args, **kwds):
        for arg in args:
            self.__dict__.update(arg)
        self.__dict__.update(kwds)
        
    def __repr__(self):
        return str(self.__dict__)

def read_header(filename):
    file = h5py.File(filename, "r")
    group = file['/Header']
    header = Struct(dict(group.attrs))
    return header
    
if __name__=="__main__":
    header = read_header("/gpfs/data/Aquila/TO/Aq-C/400/data/snapshot_127/aqc_5sig_400_ac_127.0.hdf5")
    print header.NumFilesPerSnapshot
