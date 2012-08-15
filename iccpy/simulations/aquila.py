import h5py
import glob

sim_label = { 'aqc4':'C', 'aqd4':'D', 'aqe4':'E' }

class Struct:
    def __init__(self, *args, **kwds):
        for arg in args:
            self.__dict__.update(arg)
        self.__dict__.update(kwds)
        
    def __repr__(self):
        return str(self.__dict__)

def get_dir(sim_name):
    return "/gpfs/data/Aquila/TO/Aq-%s/400/data" % sim_label[sim_name]
        
def get_data_dir(sim_name):
    return "/gpfs/data/Aquila/halo/data/Aq-%s-4" % sim_label[sim_name]
        
def get_last_snapnum(sim_name):
    if sim_name=='aqc4': return 127
    elif sim_name=='aqd4': return 511
    elif sim_name=='aqe4': return 511

def get_files(sim_name, snap_num):
    directory = get_dir(sim_name) + "/snapshot_%03d" % snapnum    
    return glob.glob(directory + "/*.hdf5")

def read_header(filename):
    file = h5py.File(filename, "r")
    group = file['/Header']
    header = Struct(dict(group.attrs))
    return header
    
def read_data(sim_name):
    file = get_data_dir(sim_name) + "/Aq-%s-4_halo_data.hdf5" % sim_label[sim_name]
    data = { name : np.array(file['/stars'][name]) for name in file['/stars'].dtype.names }
    return Struct(data)
    
def read_attr(sim_name, snap_num, part_type, attr_name):
    files = get_files(sim_name, snap_num)
    
    data = []
    for file in files:
        f = h5py.File(file, "r")
        dataset = np.array(f['/PartType%d/%s' % (part_type, attr_name)])
        
        data.append(dataset)
        
    return np.concatenate(data)
    
if __name__=="__main__":
    print read_attr("aqc4", 127, 4, 'ParticleID')

    header = read_header("/gpfs/data/Aquila/TO/Aq-C/400/data/snapshot_127/aqc_5sig_400_ac_127.0.hdf5")
    print header.NumFilesPerSnapshot
