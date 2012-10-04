import numpy as np
import os

import binary_snapshot_io

from numpy import uint32, uint64, float64, float32

_block_sizes = { "therm" : 1, "rho" : 1, "Ne" : 1, "NHI" : 1, "NHeI" : 1, "NHeIII" : 1, "sml" : 1, "pot" : 1, "accel" : 3 }

def rtype(t, swap_endian):
    if not swap_endian: 
        return t
    else:
        dt = np.dtype(t)
        return dt.newbyteorder('>')
        
def _get_files(directory, filename="", snapnum=None):
    if directory=="" and os.path.isfile(filename): return [ filename ]

    if filename=="" and snapnum==None: raise IOError("More than just directory is required to open Gadget snapshot")

    if filename!="" and snapnum==None:
        file = directory + "/" + filename
        if os.path.isfile(file): return [ file ]
        
    if filename=="" and snapnum!=None:
        #We assume that the snapshot is in a snapdir
        full_dir = directory+"/snapdir_%03d/" % snapnum        
        filelist = os.listdir(full_dir)
        
        import re
        expr = re.compile(".+_%03d.(\d+)" % snapnum)
        files = [ full_dir + name for name in filelist if expr.match(name)!=None ]
        return sorted(files, key=lambda a:int(a.split('.')[1]))
    
    #Try directory/snapshot_XXX
    file = "%s/%s_%03d" % (directory, filename, snapnum)
    if os.path.isfile(file): return [ file ]
    
    full_dir = directory+"/snapdir_%03d/" % snapnum        
    filelist = os.listdir(full_dir)
    import re
    expr = re.compile("%s_%03d.(\d+)" % (filename, snapnum))
    files = [ full_dir+name for name in filelist if expr.match(name)!=None ]
    if len(files)!=0:
        return sorted(files, key=lambda a:int(a.split('.')[1]))

    raise IOError("Unable to find Gadget snapshot file %s %s %d\n" % (directory, filename, snapnum))
             
class Header:
    def __init__(self, filename):
        """ Read the binary GADGET header file into a dictionary """
        header = binary_snapshot_io.read_snapshot_header(filename)
        self._data = header
        del self._data['buffer']
        
        self.num_particles_file = np.sum(self._data['num_particles'])
        
        if self._data['flag_doubleprecision']:
            self.dtype = rtype(float64, self._data['swap_endian'])
        else:
            self.dtype = rtype(float32, self._data['swap_endian'])
        self.dtype_width = np.dtype(self.dtype).itemsize

        #Also find length of id
        f = open(filename, mode='rb')        
        f.seek(256+16 + self.num_particles_file * 6 * self.dtype_width)
        id_block_len = np.fromfile(f, rtype(uint32, self._data['swap_endian']), 1)[0]
        if id_block_len/self.num_particles_file==4:
            self.long_ids = False
            self.id_width = 4
            self.id_type = uint32
        else:
            self.long_ids = True
            self.id_width = 8
            self.id_type = uint64            
        
        f.close()            
            
    def __str__(self):
        return self.__repr__()
            
    def __repr__(self):
        return "Gadget Snapshot Header:\n" + str(self._data)
        
    def __dir__(self):
        return sorted(set((dir(type(self)) + list(self.__dict__) + self._data.keys())))
    
    def __getattr__(self, name):
        if name not in self._data:
            raise KeyError, "'%s'" % name
        return self._data[name]
            

class Block(object):
    def __init__(self, parent, name, dim, dtype, filenames, nparts, offsets):
        self.snapshot = parent
        self.name = name
        self.dtype = dtype
        self.dtype_width = np.dtype(dtype).itemsize
        self.dim = dim
        self._loaded = np.zeros(6, dtype=np.bool)
        self._filenames = list(filenames)
        self._offsets = offsets.copy()
        self._nparts = nparts.copy()
        self._data = [ 0 ] * 6
        self._cum_nparts = np.hstack([np.zeros([len(self._filenames),1]), np.cumsum(self._nparts[:,:-1], axis=1)])
        
    def _load(self, key):
        count = 0
        #Allocate space to store data
        if self.dim==1:
            width = 1
            self._data[key] = np.empty(np.sum(self._nparts[:,key]), dtype=self.dtype)
        else:
            width = 3
            self._data[key] = np.empty([np.sum(self._nparts[:,key]), 3], dtype=self.dtype)
            
        if np.sum(self._nparts[:,key])!=0:            
            #Go through each file loading the data
            for i, filename in enumerate(self._filenames):
                #Open the file
                file = open(filename, mode='rb')
            
                #Seek to the right place in the file
                file.seek(self._offsets[i] + 4 + self.dtype_width * self._cum_nparts[i, key], 0)
            
                #Read the data from the file and place it in the data array
                data = np.fromfile(file, self.dtype, self._nparts[i, key]*width)
                
                if self.dim==1:
                    self._data[key][count:count+self._nparts[i, key]] = data
                else:
                    self._data[key][count:count+self._nparts[i, key]] = data.reshape(self._nparts[i, key], 3)
                    
                count += self._nparts[i, key]
        #print self._data[key].nbytes/1024/1024
        
        self._loaded[key] = True
      
    def __getitem__(self, key):
        if not key>=0 and key<5:
            raise IndexError, "Particle type index out of range %d" % i
    
        if not self._loaded[key]:
            self._load(key)
        return self._data[key]

class Snapshot:
    """
    A binary Gadget Snapshot
    """

    def __init__(self, directory="", filename="", snapnum=None, additional_blocks=None):
        """
        Loads a header of the snapshot
        
        Args:
            directory: the directory in which the snapshot files or directories can be found
            filename: the complete filename or the root of the snapshot filenames
            snapnum: which snapshot to load
            additional_blocks: not yet supported
        
        A directory and/or a filename must be specified. If there is more than one snapshot in the 
        directory the snapshot number must also be specified.
        """
        self._files = _get_files(directory, filename, snapnum)
        self.blocks = None
        self.header = Header(self._files[0])
        self.additional_blocks = additional_blocks
        
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        return "Gadget Snapshot"
    
    def _load_blocks(self):
        headers = [ Header(self._files[i]) for i in range(len(self._files)) ]
        self.blocks = {}
                
        nparts = np.array([ header.num_particles for header in headers ])
        nparts_per_file = np.sum(nparts, axis=1)
        
        file_offsets = np.ones(len(self._files), dtype=np.uint32) * 256 + 8
        
        self.blocks['pos'] = Block(self, "pos", 2, self.header.dtype, self._files, nparts, file_offsets)
        file_offsets += 3 * self.header.dtype_width * nparts_per_file + 8
        
        self.blocks['vel'] = Block(self, "vel", 2, self.header.dtype, self._files, nparts, file_offsets)
        file_offsets += 3 * self.header.dtype_width * nparts_per_file + 8
        
        self.blocks['id'] = Block(self, "id", 1, self.header.id_type, self._files, nparts, file_offsets)
        file_offsets += self.header.id_width * nparts_per_file + 8

        #Mass blocks are evil
        masses = np.array([ header.mass for header in headers ])
        nparts_mass = np.zeros_like(nparts)
        nparts_mass[np.where(masses==0)] = nparts[np.where(masses==0)]
        
        if np.sum(nparts_mass)!=0:
            self.blocks['mass'] = Block(self, "mass", 1, self.header.dtype, self._files, nparts_mass, file_offsets)
            file_offsets += self.header.dtype_width * np.sum(nparts_mass, axis=1) + 8
            #Need to set non-read masses
            for i in range(6):
                if np.sum(masses[:,i])!=0: 
                    self.blocks['mass']._loaded[i] = True
                    self.blocks['mass']._data[i] = np.array(masses[np.where(masses[:,i]!=0)[0][0],i])
        else:
            self.blocks['mass'] = self.header.mass
            
        for name in self.additional_blocks:
            if name not in _block_sizes:
                raise KeyError, "Unknown block name %s" % name
            
            dim = 1 if _block_sizes[name]==1 else 2
            self.blocks[name] = Block(self, name, dim, self.header.id_type, self._files, nparts, file_offsets)
            file_offsets +=  _block_sizes[name] * self.header.id_width * nparts_per_file + 8
        
    def __getattr__(self, name):        
        if self.blocks is None:
            self._load_blocks()
    
        if name not in self.blocks:
            raise KeyError, "Unknown block name %s" % name
        
        return self.blocks[name]
        
def snapshot_files(directory, filename="", snapnum=None, additional_blocks=None):
    files = _get_files(directory, filename, snapnum)

    for file in files:
        yield Snapshot("", file)

def save(filename, header, pos, vel, ids=None, masses=None, extra_data=None):
    if ids is None:
        ids = np.arange(pos.shape[0], dtype=np.uint32)
    
    if header is None:
        if masses is None:
            mass_header = np.array([0, 1, 0, 0, 0, 0])
        else:
            mass_header = np.array([0, 0, 0, 0, 0, 0])
    
        num_particles = np.array([0, pos.shape[0], 0, 0, 0, 0])
        header = dict((('num_particles', num_particles),
                   ('mass', mass_header), ('time',0.0), ('redshift',0), 
                   ('flag_sfr',0) , ('flag_feedback',0), 
                   ('num_particles_total', num_particles), 
                   ('flag_cooling',0), ('num_files',1), ('boxsize',0.0), \
                   ('omega0', 0.0), ('omegaLambda', 0.0), ('hubble0', 0.0), ('flag_stellarage',0), \
                   ('buffer', [0]*56), ('flag_metals', 0), ('npartTotalHighWord', [0,0,0,0,0,0]), 
                   ('flag_entropy_instead_u', 0), ('flag_doubleprecision', 1)))
                   
        binary_snapshot_io.write_snapshot_file(filename, header, pos, vel, ids, masses, extra_data)
    elif isinstance(header, Header):
        header._data['buffer'] = np.empty(56, np.uint8)
        binary_snapshot_io.write_snapshot_file(filename, header._data, pos, vel, ids, masses, extra_data)    
        del header._data['buffer']
    else:
        binary_snapshot_io.write_snapshot_file(filename, header, pos, vel, ids, masses, extra_data)
