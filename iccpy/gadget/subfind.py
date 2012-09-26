import numpy as np
import time
import iccpy.utils
import os.path

class Group:
    def __init__(self, parent, idx, block):
        self._parent = parent
        self._idx = idx
        self._block = block
        self.num_particles = block["num_particles"][idx]
        self.offset = block["offset"][idx]        
        
    def __getattr__(self, name):
        if name=="id":
            #Only load the ids if we need them
            if self._parent.ids is None: self._parent._load_ids()
            #Return a slice of the parent array
            return self._parent.ids[self.offset:self.offset+self.num_particles]
        else if name in self._block:
            return self._block[name][self._idx]
        else:
            raise KeyError, "Unknown attribute name %s" % name
    
def find_group_membership(ids, groupids, length, offset, find_mbrank=False):
    """
    Given a subfind output determine group membership for each particle
    in the snapshot. Returns -1 for particles not in a group.
    """

    # Get unique version of the ids to find
    (uids, revidx) = np.unique(ids, return_inverse=True)

    # Find matching elements between ids and groupids
    ptr = iccpy.utils.match(groupids, uids)

    # Find group membership for each element in groupids
    grnr   = np.zeros(len(groupids), dtype=int) - 1
    mbrank = np.zeros(len(groupids), dtype=int)
    for igroup in range(len(length)):
        grnr[offset[igroup]:offset[igroup]+length[igroup]]   = igroup
        mbrank[offset[igroup]:offset[igroup]+length[igroup]] = arange(length[igroup], dtype=int) + 1

    # Find group membership for each particle in the snapshot
    snap_grnr   = np.zeros(len(uids), dtype=int) - 1
    snap_mbrank = np.zeros(len(uids), dtype=int)
    ind = (ptr >= 0)
    snap_grnr[ptr[ind]]   = grnr[ind]
    snap_mbrank[ptr[ind]] = mbrank[ind]

    # Return reconstructed arrays including results for duplicate IDs
    if find_mbrank:
        return (snap_grnr[revidx], snap_mbrank[revidx])
    else:
        return snap_grnr[revidx]

class SubfindCatalogue:
    def __init__(self, directory, snapshot_num, float_type=np.float32,
                 SO_VEL_DISPERSIONS=False, SO_BAR_INFO=False):
        self.directory = directory
        self.snapshot_num = snapshot_num
        self.fof_group = None
        self.subhalo = None

        self.ids = None
        self.sub_lengths = None
        self.sub_offsets = None
        self.grp_lengths = None
        self.grp_offsets = None

        self.id_size = self._determine_id_size()
        
        num_files = self._determine_num_files()
        print num_files
        
    def _determine_num_files(self):
        filename = "%s/groups_%03d/subhalo_tab_%03d.0" % (self.directory, self.snapshot_num, self.snapshot_num)
        if not os.path.exists(filename): raise IOError("Unable to find subtab file %s" % filename)

        f = open(filename, mode='rb')
        f.seek(20, 1)        
        num_files = np.fromfile(f, u32, 1)[0]
        f.close()
        
        byteswap = (num_files<0) or (num_files>=65536)
        if byteswap:
            num_files.byteswap()

        return num_files
                
    def _determine_id_size(self):
        filename = "%s/groups_%03d/subhalo_ids_%03d.0" % (self.directory, self.snapshot_num, self.snapshot_num)
        if not os.path.exists(filename): return 4
        
        f = open(filename, mode='rb')
    
        f.seek(8, 1)
        num_ids       = np.fromfile(f, np.uint32, 1)[0]
        num_ids_total = np.fromfile(f, np.uint64, 1)[0]
        num_files     = np.fromfile(f, np.uint32, 1)[0]
        
         # Assume large/negative ntask means wrong endian!
        byteswap = (num_files<0) or (num_files>=65536)
        
        if byteswap: num_ids.byteswap()    
        f.close()

        filelength = os.path.getsize(filename)
    
        if filelength == 28 + 8 * num_ids: return 8
        elif filelength == 28 + 4 * num_ids: return 4
        else: raise Exception('Unable to determine size of ID type from file length')
        
    def _load_ids(self):
        filebase = "%s/groups_%03d/subhalo_ids_%03d" % (self.directory, self.snapshot_num, self.snapshot_num)
        
        filename = "%s.0" % filebase
        if not os.path.exists(filename): raise IOError("Unable to find ID files %s" % filebase)
        f = open(filename, mode='rb')
        
        #Just read the information we need
        f.seek(8, 1)
        num_ids       = np.fromfile(f, np.uint32, 1)[0]
        num_ids_total = np.fromfile(f, np.uint64, 1)[0]
        num_files     = np.fromfile(f, np.uint32, 1)[0]
        
         # Assume large/negative ntask means wrong endian!
        byteswap = (num_files<0) or (num_files>=65536)
        
        if byteswap: 
            num_ids.byteswap()
            num_files.byteswap()
        f.close()

        id_list = []
        id_type = np.uint32 if self.id_size==4: else np.uint64
    
        for i in range(num_files):
            filename = "%s.%d" % (filebase, i)
            f = open(filename, mode='rb')
        
            f.seek(8, 1)
            num_ids = np.fromfile(f, np.uint32, 1)[0]
            if byteswap: num_ids.byteswap()
            
            f.seek(16, 1)
            ids = np.fromfile(f, id_type, num_ids)
            f.close()
    
            id_list.append(ids)
    
        self.ids = np.concatenate(id_list)
        if byteswap: self.ids.byteswap()
        
    def _read_subtab_file(self, fname, float_type=np.float32, id_size=4, SO_VEL_DISPERSIONS=False,
                          SO_BAR_INFO=False):
        f = open(fname,"r")
        # Header
        ngroups    = np.fromfile(f, dtype=np.int32, count=1)[0]
        totngroups = np.fromfile(f, dtype=np.int32, count=1)[0]
        nids       = np.fromfile(f, dtype=np.int32, count=1)[0]
        totnids    = np.fromfile(f, dtype=np.int64, count=1)[0]
        ntask      = np.fromfile(f, dtype=np.int32, count=1)[0]
        nsubgroups = np.fromfile(f, dtype=np.int32, count=1)[0]
        totnsubgroups = np.fromfile(f, dtype=np.int32, count=1)[0]

        # Assume large/negative ntask means wrong endian!
        byteswap = (ntask<0) or (ntask>=65536)

        # Byteswap header if necessary
        if byteswap:
            ngroups    = ngroups.byteswap()
            totngroups = totngroups.byteswap()
            nids       = nids.byteswap()
            totnids    = totnids.byteswap()
            ntask      = ntask.byteswap()
            nsubgroups    = nsubgroups.byteswap()
            totnsubgroups = totnsubgroups.byteswap()

        fof_block = {}
        subhalo_block = {}
        
        # FoF group properties
        fof_block["num_particles"]   = np.fromfile(f, dtype=np.int32,   count=ngroups)
        fof_block["offset"]          = np.fromfile(f, dtype=np.int32,   count=ngroups)
        fof_block["mass"]            = np.fromfile(f, dtype=float_type, count=ngroups)
        fof_block["pot_min"]         = np.fromfile(f, dtype=float_type, count=3*ngroups).reshape((ngroups,3))
                
        fof_block["mass_mean200"]    = np.fromfile(f, dtype=float_type, count=ngroups)
        fof_block["radius_mean200"]  = np.fromfile(f, dtype=float_type, count=ngroups)
        fof_block["mass_crit200"]    = np.fromfile(f, dtype=float_type, count=ngroups)
        fof_block["radius_crit200"]  = np.fromfile(f, dtype=float_type, count=ngroups)
        fof_block["mass_tophat"]     = np.fromfile(f, dtype=float_type, count=ngroups)
        fof_block["radius_tophat_"]  = np.fromfile(f, dtype=float_type, count=ngroups)

        if SO_VEL_DISPERSIONS:
            fof_block["veldisp_mean200"] = np.fromfile(f, dtype=float_type, count=ngroups)
            fof_block["veldisp_crit200"] = np.fromfile(f, dtype=float_type, count=ngroups)
            fof_block["veldisp_tophat"]  = np.fromfile(f, dtype=float_type, count=ngroups)

        if SO_BAR_INFO:
            fof_block["gas_mass"]  = np.fromfile(f, dtype=float_type, count=3*ngroups).reshape((ngroups,3))
            fof_block["star_mass"] = np.fromfile(f, dtype=float_type, count=3*ngroups).reshape((ngroups,3))
            fof_block["gas_temp"]  = np.fromfile(f, dtype=float_type, count=3*ngroups).reshape((ngroups,3))
            fof_block["gas_xlum"]  = np.fromfile(f, dtype=float_type, count=3*ngroups).reshape((ngroups,3))

        fof_block["contaimination_count"] = np.fromfile(f, dtype=np.int32,   count=ngroups)
        fof_block["contaimination_mass"]  = np.fromfile(f, dtype=float_type, count=ngroups)
        fof_block["num_subhaloes"]        = np.fromfile(f, dtype=np.int32,   count=ngroups)
        fof_block["first_subhalo_offset"] = np.fromfile(f, dtype=np.int32,   count=ngroups)

        # Subhalo properties
        subhalo_block["num_particles"]       = np.fromfile(f, dtype=np.int32,   count=nsubgroups)
        subhalo_block["offset"]       = np.fromfile(f, dtype=np.int32,   count=nsubgroups)
        subhalo_block["parent"]       = np.fromfile(f, dtype=np.int32,   count=nsubgroups)
        subhalo_block["mass"]      = np.fromfile(f, dtype=float_type, count=nsubgroups)
        subhalo_block["pot_min"]       = np.fromfile(f, dtype=float_type, count=3*nsubgroups).reshape((nsubgroups,3))
        subhalo_block["velocity"]       = np.fromfile(f, dtype=float_type, count=3*nsubgroups).reshape((nsubgroups,3))
        subhalo_block["com"]      = np.fromfile(f, dtype=float_type, count=3*nsubgroups).reshape((nsubgroups,3))
        subhalo_block["spin"]         = np.fromfile(f, dtype=float_type, count=3*nsubgroups).reshape((nsubgroups,3))
        subhalo_block["vel_disp"]      = np.fromfile(f, dtype=float_type, count=nsubgroups)
        subhalo_block["vel_max"]         = np.fromfile(f, dtype=float_type, count=nsubgroups)
        subhalo_block["radius_vel_max"]        = np.fromfile(f, dtype=float_type, count=nsubgroups)
        subhalo_block["half_mass_radius"]        = np.fromfile(f, dtype=float_type, count=nsubgroups)

        if id_size == 4:
            subhalo_block["most_bound_particle_id"] = np.fromfile(f, dtype=np.int32, count=nsubgroups)
        elif id_size == 8:
            subhalo_block["most_bound_particle_id"] = np.fromfile(f, dtype=np.int64, count=nsubgroups)
        else:
            raise Exception("id_size must be 4 or 8!")

        subhalo_block["grnr"] = np.fromfile(f, dtype=np.int32, count=nsubgroups)

        # Byteswap data if necessary
        if byteswap:
            for key in fof_block.keys():
                fof_block[key] = fof_block[key].byteswap()
            for key in subhalo_block.keys():
                subhalo_block[key] = subhalo_block[key].byteswap()                
                
        #Convert into set of fofgroups and subhaloes
        groups = [ Group(self, i, fof_block) for i in range(ngroups) ]
        subhaloes = [ Group(self, i, subhalo_block) for i in range(ngroups) ]
        
        return groups, subhaloes
    
     def subhalo_index(self, ids, find_mbrank=False):
        """Calculate subhalo membership for the given particle IDs"""
        if self.ids is None: self._load_ids()
        res = find_group_membership(ids, self.ids, self.sub_lengths, self.sub_offsets, find_mbrank)
        return res

    def fofgroup_index(self, ids):
        """Calculate FoF group membership for the given particle IDs"""
        if self.ids is None: self._load_ids()
        res = find_group_membership(ids, self.ids, self.grp_lengths, self.grp_offsets)
        return res
