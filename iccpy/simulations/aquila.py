import h5py
import glob
import numpy as np

sim_label = { 'aqc4':'C', 'aqd4':'D', 'aqe4':'E' }
max_id = { 'aqc4' : (53967599, 103), 'aqd4' : (41535584, 61), 'aqe4' : (9223281412051116774, 80) }
last_snapnum = { 'aqc4' : 127, 'aqd4' : 511, 'aqe4' : 511 }

GAS = 0
DARK_MATTER = DM = 1
STARS = 4

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

def get_fig_dir(sim_name):
    return "/gpfs/data/Aquila/halo/figs/Aq-%s-4" % sim_label[sim_name]

def get_files(sim_name, snap_num):
    directory = get_dir(sim_name) + "/snapshot_%03d" % snap_num    
    return glob.glob(directory + "/*.hdf5")

def read_header(filename):
    file = h5py.File(filename, "r")
    group = file['/Header']
    header = Struct(dict(group.attrs))
    return header
    
def read_data(sim_name):
    filename = get_data_dir(sim_name) + "/Aq-%s-4_halo_data.hdf5" % sim_label[sim_name]
    file = h5py.File(filename, 'r')
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
    
def make_unique_star_ids(sim_name, snap_num, ids):
    child_ids = read_attr(sim_name, snap_num, 4,'ChildIDforStars')
    unique_ids = max_id[sim_name][0] + (max_id[sim_name][1] + 1) * (1 + ids) + child_ids
    return unique_ids
    
def get_halo_centre(sim_name, snap_num):
    fname = get_dir(sim_name) + "/groups_%03d/subhalo_tab_%03d.0" % (snap_num, snap_num)
    file = open(fname, "rb")
    
    #Load number of groups and subgroups
    num_groups = np.fromfile(file, np.int32, 1)[0]
    file.seek(20, 1)
    num_subgroups = np.fromfile(file, np.int32, 1)[0]
    file.seek(4 + 76*num_groups + 16*num_subgroups, 1)
    
    centre = np.fromfile(file, np.float32, 3)
    file.close()
    return centre

if __name__=="__main__":
    sim_name = "aqd4"
    #make_unique_star_ids(sim_name, last_snapnum[sim_name], 

    #ids = read_attr("aqd4", 127, 4, 'ParticleIDs')
    #print np.where(ids==588263426)

    #data = read_data("aqd4")
    #idxs = np.where((data.origin==0) & (data.snap_last_in_sub==127))
    #ids_snapshot = data.ID[idxs]

    #print ids_snapshot
