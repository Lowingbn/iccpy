import binary_group_io
import snapshot
import itertools
import numpy as np

class FoFGroup:
    """Store information about each FoF group"""
    def __init__(self, snapshot_num):
        self.snapshot_num = snapshot_num
        self.npart = 0
        self.ids = None
    

class SubfindGroup:
    """Store information about each subfind group"""
    def __init__(self, snapshot_num):
        self.snapshot_num = snapshot_num
        self.npart = 0
        self.pot_min = None
        self.most_bound_particle = None
        self.mass = 0
        self.half_mass_radius = 0
        self.v_max = 0
        self.ids = None
        
def load_FoF_groups(directory, snapshot_num, loadIDs=False, endianness='native'):
    if loadIDs:
        filebase = "%s/groups_%03d/group_ids_%03d" % (directory, snapshot_num, snapshot_num)
        ids = binary_group_io.read_IDs(filebase, endianness)
    else:
        ids = None
        
    #Load first file
    filename = "%s/groups_%03d/group_tab_%03d.0" % (directory, snapshot_num, snapshot_num)
    props = binary_group_io.read_group_file(filename, None, endianness)[0]
    
    groups = []
    
    for i in range(props['num_files']):
        filename = "%s/groups_%03d/group_tab_%03d.%d" % (directory, snapshot_num, snapshot_num, i)
        groups_i = binary_group_io.read_group_file(filename, ids, endianness)[1]
        
        groups.extend(convert_to_fof_groups(groups_i, snapshot_num))
        
    return groups
            
def load_subfind_groups(directory, snapshot_num, loadIDs=False, endianness='native'):
    if loadIDs:
        filebase = "%s/groups_%03d/subhalo_ids_%03d" % (directory, snapshot_num, snapshot_num)
        ids = binary_group_io.read_IDs(filebase, endianness)
    else:
        ids = None
        
    #Load first file
    filename = "%s/groups_%03d/subhalo_tab_%03d.0" % (directory, snapshot_num, snapshot_num)
    props = binary_group_io.read_subfind_file(filename, None, endianness)[0]
    
    groups = []
    subgroups = []
    
    for i in range(props['num_files']):
        filename = "%s/groups_%03d/subhalo_tab_%03d.%d" % (directory, snapshot_num, snapshot_num, i)
        groups_i, subgroups_i = binary_group_io.read_subfind_file(filename, ids, endianness)[1:3]
        
        groups.extend(convert_to_fof_groups(groups_i, snapshot_num))
        subgroups.extend(convert_to_subfind_groups(subgroups_i, snapshot_num))
        
    return groups, subgroups

def convert_to_fof_groups(groups_props, snapshot_num):
    """Convert a dictionary of group properties into a list of group objects"""
    groups = []
    
    for i in range(len(groups_props['npart'])):
        group = FoFGroup(snapshot_num)
        group.npart = groups_props['npart'][i]
        
        if 'ids' in groups_props:
            group.ids = groups_props['ids'][i]
        
        groups.append(group)
    
    return groups

def convert_to_subfind_groups(subgroups_props, snapshot_num):
    subgroups = []
    
    for i in range(len(subgroups_props['npart'])):
        subgroup = SubfindGroup(snapshot_num)
        subgroup.npart = subgroups_props['npart'][i]
        subgroup.pot_min = subgroups_props['pot_min'][i]
        subgroup.most_bound_particle = subgroups_props['most_bound_particle'][i]
        subgroup.mass = subgroups_props['mass'][i]
        subgroup.v_max = subgroups_props['vel_max'][i]
        subgroup.half_mass_radius = subgroups_props['half_mass_radius'][i]
        subgroup.vel = subgroups_props['vel'][i]
        
        if 'ids' in subgroups_props:
            subgroup.ids = subgroups_props['ids'][i]
        
        subgroups.append(subgroup)
    
    return subgroups

def get_subgroup_idx(id, directory, snapshot_num, ids=None, endianness='native'):
    if ids is None:
        filebase = "%s/groups_%03d/subhalo_ids_%03d" % (directory, snapshot_num, snapshot_num)
        ids = binary_group_io.read_IDs(filebase, endianness)
        
        idx = np.flatnonzero(ids==id)[0]
        del ids
    else:
        idx = np.flatnonzero(ids==id)[0]
    
    #Open file until we find the right one
    subgroup_counter = 0
    for i in itertools.count():
        filename = "%s/groups_%03d/subhalo_tab_%03d.%d" % (directory, snapshot_num, snapshot_num, i)
        props, groups, subgroups = binary_group_io.read_subfind_file(filename, endianness)
        
        for j in range(len(subgroups['npart'])):
            if subgroups['offsets'][j]<=idx and subgroups['offsets'][j]+subgroups['npart'][j]>idx:
                del groups
                del subgroups
                
                return subgroup_counter + j
        
        subgroup_counter += props['num_subgroups']

        del groups
        del subgroups

def get_subgroup_ids(subgroup_num, directory, snapshot_num, ids=None, endianness='native'):
    if ids is None:
        filebase = "%s/groups_%03d/subhalo_ids_%03d" % (directory, snapshot_num, snapshot_num)
        ids = binary_group_io.read_IDs(filebase, endianness)
        
    try:
        iterator = iter(subgroup_num)
        selectedIDs = set()
        for subgroup_num in iterator:
            selectedIDs.update(get_subgroup_ids(subgroup_num, directory, snapshot_num, ids))
        return selectedIDs
    except TypeError: 
        #Open file until we find the right one
        for i in itertools.count():
            filename = "%s/groups_%03d/subhalo_tab_%03d.%d" % (directory, snapshot_num, snapshot_num, i)
            props, groups, subgroups = binary_group_io.read_subfind_file(filename, endianness)
        
            if subgroup_num<props['num_subgroups']: break
        
            subgroup_num -= props['num_subgroups']
            del groups
            del subgroups
    
        offset = subgroups['offsets'][subgroup_num]
        len = subgroups['npart'][subgroup_num]
    
        del groups
        del subgroups
    
        return set(ids[offset:offset+len])
    
def get_particles(directory, filename, snapshot_num, select_ids):
    ids, poses, vels = [], [], []
    
    for header, res in snapshot.load_snapshot_files(directory, filename, snapshot_num):
        mass = res['mass'][1]
        for id, pos, vel in itertools.izip(res['id'], res['pos'], res['vel']):
            if id in select_ids:
                ids.append(id)
                poses.append(pos)
                vels.append(vel)
    
        del res
    
    return np.array(ids), np.array(poses), np.array(vels), np.array([mass] * len(ids))

def get_subgroup_particles(directory, filename, snapshot_num, subgroup_nums):
    selectedIDs = get_subgroup_ids(subgroup_nums, directory, snapshot_num)
    
    return get_particles(directory, filename, snapshot_num, selectedIDs)
