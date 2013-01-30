import numpy as np
import os
import h5py
import re

import binary_snapshot_io

from numpy import uint32, uint64, float64, float32

_block_element_nums = { "therm" : 1, "rho" : 1, "Ne" : 1, "NHI" : 1, "NHeI" : 1, "NHeIII" : 1, "sml" : 1, "pot" : 1, "accel" : 3 }

def _rtype(t, swap_endian):
    if not swap_endian: 
        return t
    else:
        dt = np.dtype(t)
        return dt.newbyteorder('>')
                
def _convert(name):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()        
        
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
        
        expr = re.compile(".+_%03d.(\d+)" % snapnum)
        files = [ full_dir + name for name in filelist if expr.match(name)!=None ]
        return sorted(files, key=lambda a:int(a.split('.')[1]))
    
    #Try directory/snapshot_XXX
    if directory!="" and filename!="" and snapnum!=None:
        file = "%s/%s_%03d" % (directory, filename, snapnum)
        if os.path.isfile(file): return [ file ]
    
        full_dir = directory+"/snapdir_%03d/" % snapnum        
        filelist = os.listdir(full_dir)
        expr = re.compile("%s_%03d.(\d+)" % (filename, snapnum))
        files = [ full_dir+name for name in filelist if expr.match(name)!=None ]
        if len(files)!=0:
            return sorted(files, key=lambda a:int(a.split('.')[1]))

    if snapnum==None:
        errorString = "%s %s" % (directory, filename)
    else:
        errorString = "%s %s %d" % (directory, filename, snapnum)        
    raise IOError("Unable to find Gadget snapshot file %s\n" % errorString)
             
class GadgetBinaryHeader:
    def __init__(self, filename, format=1):
        """ Read the binary GADGET header """
        file = open(filename, mode='rb')
        if format==2: 
            file.seek(16) # skip 4 HEAD 264 4
        header = binary_snapshot_io.read_header(file)
        file.close()
        
        self._data = header
        del self._data['buffer']
        
        self.num_particles_file = np.sum(self._data['num_particles'])
        
        if self._data['flag_doubleprecision']:
            self.dtype = _rtype(float64, self._data['swap_endian'])
        else:
            self.dtype = _rtype(float32, self._data['swap_endian'])
        self.dtype_width = np.dtype(self.dtype).itemsize

        #Also find length of id
        f = open(filename, mode='rb')        
        f.seek(256+24 + int(self.num_particles_file) * 6 * int(self.dtype_width))
        if format==2:
            f.seek(16*4,1) # skip marker prefixes (HEAD, POS, VEL, IDS)

        id_block_len = np.fromfile(f, _rtype(uint32, self._data['swap_endian']), 1)[0]
        if id_block_len/self.num_particles_file==4:
            self.long_ids = False
            self.id_width = 4
            self.id_type = _rtype(uint32, self._data['swap_endian'])
        else:
            self.long_ids = True
            self.id_width = 8
            self.id_type = _rtype(uint64, self._data['swap_endian'])            
        
        f.close()            
            
    def __str__(self):
        return self.__repr__()
            
    def __repr__(self):
        return "Gadget Binary Snapshot Header:\n" + str(self._data)
        
    def __dir__(self):
        return sorted(set((dir(type(self)) + list(self.__dict__) + self._data.keys())))
    
    def __getattr__(self, name):
        if name not in self._data:
            raise KeyError, "'%s'" % name
        return self._data[name]
            

class GadgetBinaryBlock(object):
    def __init__(self, parent, name, elements, dtype, filenames, nparts, offsets):
        self.snapshot = parent
        self.name = name
        self.dtype = dtype
        self.dtype_width = np.dtype(dtype).itemsize
        self.elements = int(elements)
        self._loaded = np.zeros(6, dtype=np.bool)
        self._filenames = list(filenames)
        self._offsets = offsets.copy()
        self._nparts = nparts.copy()
        self._data = [ 0 ] * 6
        self._cum_nparts = np.hstack([np.zeros([len(self._filenames),1], dtype=self._nparts.dtype), np.cumsum(self._nparts[:,:-1], axis=1, dtype=np.int64)])
        
    def _load(self, key):

        #Allocate space to store data
        if self.elements==1:
            self._data[key] = np.empty(np.sum(self._nparts[:,key]), dtype=self.dtype)
        else:
            self._data[key] = np.empty([np.sum(self._nparts[:,key]), self.elements], dtype=self.dtype)
            
        if np.sum(self._nparts[:,key])!=0:            
            #Go through each file loading the data
            count = 0
            for i, filename in enumerate(self._filenames):
                #Open the file
                file = open(filename, mode='rb')
            
                #Seek to the right place in the file
                file.seek(self._offsets[i] + 4 + self.dtype_width * self._cum_nparts[i, key] * self.elements, 0)
            
                #Read the data from the file and place it in the data array
                data = np.fromfile(file, self.dtype, self._nparts[i, key]*self.elements)

                if self.elements==1:
                    self._data[key][count:count+self._nparts[i, key]] = data
                else:
                    self._data[key][count:count+self._nparts[i, key]] = data.reshape(self._nparts[i, key],
                                                                                     self.elements)
                file.close()

                count += self._nparts[i, key]
        #print self._data[key].nbytes/1024/1024
        
        self._loaded[key] = True

      
    def __getitem__(self, key):
        if not key>=0 and key<5:
            raise IndexError, "Particle type index out of range %d" % i
    
        if not self._loaded[key]:
            self._load(key)
        return self._data[key]
        
##########################################################################################

def _determine_binary_format(filename):
    f = open(filename, mode='rb')
    r = np.fromfile(f, np.uint32, 1)[0]
    f.close()
    if r==8 or r==134217728:
        return 2
    elif r==256 or r==65536:
        return 1
    else:
        raise IOError, "File corrupt. First integer is: " + str(r)

class GadgetBinaryFormat1Snapshot:
    """
    A binary format 1 Gadget Snapshot
    """

    def __init__(self, directory="", filename="", snapnum=None, files=None, additional_blocks=None):
        """
        Creates a Gadget binary snapshot object and loads the header
        
        Keyword arguments:
            directory          -- the directory in which the snapshot files or directories can be found
            filename:          -- the complete filename or the root of the snapshot filenames
            snapnum:           -- which snapshot to load
            files:             -- a list of files to load
            additional_blocks  -- a list of blocks other than the basic that exist in the snapshot
        
        A directory and/or a filename or a list of files must be specified. If there is more than one snapshot in the 
        directory the snapshot number must also be specified.
        """
        if files is None:
            self._files = _get_files(directory, filename, snapnum)
        else:
            self._files = list(files)
        self._blocks = None
        self.header = GadgetBinaryHeader(self._files[0])
        self._additional_blocks = additional_blocks
        
        if len(files)==1:
            self.num_particles = self.header.num_particles
        else:
            self.num_particles = self.header.num_particles_total
        
        del self.header._data['num_particles']
        del self.header.num_particles_file
        
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        return "Gadget Binary Format 1 Snapshot"
        
    def __dir__(self): 
        if self._blocks is None:
            self._load_blocks()   
        return sorted(set((list(self.__dict__) + self._blocks.keys())))         
    
    def _load_blocks(self):
        headers = [ GadgetBinaryHeader(self._files[i]) for i in range(len(self._files)) ]
        self._blocks = {}
                
        nparts = np.array([ header.num_particles for header in headers ], dtype=np.int64)
        nparts_per_file = np.sum(nparts, axis=1)
        
        file_offsets = np.ones(len(self._files), dtype=np.uint32) * 256 + 8
        
        self._blocks['pos'] = GadgetBinaryBlock(self, "pos", 3, self.header.dtype, self._files, nparts, file_offsets)
        file_offsets += 3 * self.header.dtype_width * nparts_per_file + 8
        
        self._blocks['vel'] = GadgetBinaryBlock(self, "vel", 3, self.header.dtype, self._files, nparts, file_offsets)
        file_offsets += 3 * self.header.dtype_width * nparts_per_file + 8
        
        self._blocks['id'] = GadgetBinaryBlock(self, "id", 1, self.header.id_type, self._files, nparts, file_offsets)
        file_offsets += self.header.id_width * nparts_per_file + 8

        #Mass blocks are evil
        masses = np.array([ header.mass for header in headers ])
        nparts_mass = np.zeros_like(nparts)
        nparts_mass[np.where(masses==0)] = nparts[np.where(masses==0)]

        if np.sum(nparts_mass)!=0:
            self._blocks['mass'] = GadgetBinaryBlock(self, "mass", 1, self.header.dtype, self._files, nparts_mass, file_offsets)
            file_offsets += self.header.dtype_width * np.sum(nparts_mass, axis=1) + 8
            #Need to set non-read masses
            for i in range(6):
                if np.sum(masses[:,i])!=0: 
                    self._blocks['mass']._loaded[i] = True
                    self._blocks['mass']._data[i] = np.ones(headers[0].num_particles_total[i]) * headers[0].mass[i]  #np.array(masses[np.where(masses[:,i]!=0)[0][0],i])
        else:
            self._blocks['mass'] = self.header.mass
        
        if self._additional_blocks is not None:
            for name in self._additional_blocks:
                if name not in _block_element_nums:
                    raise KeyError, "Unknown block name %s" % name
            
                dim = 1 if _block_element_nums[name]==1 else 2
                self._blocks[name] = GadgetBinaryBlock(self, name, dim, self.header.dtype, self._files, nparts, file_offsets)
                file_offsets +=  _block_element_nums[name] * self.header.dtype_width * nparts_per_file + 8
        
    def __getattr__(self, name):        
        if self._blocks is None:
            self._load_blocks()
    
        if name not in self._blocks:
            raise KeyError, "Unknown block name %s" % name
        
        return self._blocks[name]
        
class GadgetBinaryFormat2Snapshot:
    """
    A binary format 2 Gadget Snapshot
    """

    def __init__(self, directory="", filename="", snapnum=None, files=None, label_table=None):
        """
        Creates a Gadget binary snapshot object and loads the header
        
        Keyword arguments:
            directory          -- the directory in which the snapshot files or directories can be found
            filename:          -- the complete filename or the root of the snapshot filenames
            snapnum:           -- which snapshot to load
            files:             -- a list of files to load
            label_table:      -- a dictionary of additional labels (e.g. {"ACCE":(3, (1,1,1,1,1,1))})
                                  which specifies the number of fields and which particles this label is for.

        A directory and/or a filename or a list of files must be specified. If there is more than one snapshot in the 
        directory the snapshot number must also be specified.
        """
        if files is None:
            self._files = _get_files(directory, filename, snapnum)
        else:
            self._files = list(files)
        self._blocks = None
        self.header = GadgetBinaryHeader(self._files[0], format=2)
        
        if len(files)==1:
            self.num_particles = self.header.num_particles
        else:
            self.num_particles = self.header.num_particles_total
        
        del self.header._data['num_particles']
        del self.header.num_particles_file
        
        self._label_table = label_table
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        return "Gadget Binary Format 2 Snapshot"
        
    def __dir__(self): 
        if self._blocks is None:
            self._load_blocks()
        return sorted(set((list(self.__dict__) + self._blocks.keys())))
        
    def _load_blocks(self):
        headers = [ GadgetBinaryHeader(self._files[i], format=2) for i in range(len(self._files)) ]
        block_keys = _get_all_format2_names(self._files[0], headers[0].swap_endian)
        self._blocks = {}
        print 'Keys in file', block_keys
        nparts = np.array([ header.num_particles for header in headers ])
        nparts_per_file = np.sum(nparts, axis=1)
        
        file_offsets = np.ones(len(self._files), dtype=np.uint32) * 256 + 8 + 32
        
        self._blocks[block_keys[1]] = GadgetBinaryBlock(self, "pos", 3, self.header.dtype, self._files, nparts, file_offsets)
        file_offsets += 3 * self.header.dtype_width * nparts_per_file + 8 + 16
        
        self._blocks[block_keys[2]] = GadgetBinaryBlock(self, "vel", 3, self.header.dtype, self._files, nparts, file_offsets)
        file_offsets += 3 * self.header.dtype_width * nparts_per_file + 8 + 16
        
        self._blocks[block_keys[3]] = GadgetBinaryBlock(self, "id", 1, self.header.id_type, self._files, nparts, file_offsets)
        file_offsets += self.header.id_width * nparts_per_file + 8 + 16
        
        #Mass blocks are evil
        masses = np.array([ header.mass for header in headers ])
        nparts_mass = np.zeros_like(nparts)
        nparts_mass[np.where(masses==0)] = nparts[np.where(masses==0)]

        if np.sum(nparts_mass)!=0:
            assert(block_keys[4]=='MASS')  # should be mass
            mass = block_keys[4]

            self._blocks[mass] = GadgetBinaryBlock(self, "mass", 1, self.header.dtype, self._files, nparts_mass, file_offsets)
            file_offsets += self.header.dtype_width * np.sum(nparts_mass, axis=1) + 8 + 16
            #Need to set non-read masses
            for i in range(6):
                if np.sum(masses[:,i])!=0:
                    self._blocks[mass]._loaded[i] = True
                    self._blocks[mass]._data[i] = np.ones(headers[0].num_particles_total[i]) * headers[0].mass[i]  
        else:
            # there is no mass block
            self._blocks[mass] = self.header.mass
        
        for i in range(5, len(block_keys)):
            name = block_keys[i]
            if name not in self._label_table:
                raise KeyError, "<%s> is not a known label (try adding to the label_table?)"%name 

            elements, used_parts = self._label_table[name]

            ## TODO this next assertion should probably go in to a "GadgetBinaryBlockFormat2"
            ## to be cleaner.
            f = open(self._files[0], 'rb')
            f.seek(file_offsets[0]-12)
            found_name = f.read(4)
            #num_bytes = np.fromfile(f, uint32, 1)
            #print 'bytes to read', num_bytes
            f.close()
            if name!=found_name:
                raise Exception("Expected <%s> in file but found <%s> (perhaps due to incorrect expected length of records for <%s> ?"%(name, found_name, block_keys[i-1]))
            f.close()
            #######



            nparts_var = nparts * used_parts #[npt * present for npt, present in zip(nparts, (1,0,0,0,0,0))]
            nparts_var_per_file = np.sum(nparts_var, axis=1)
            #print 'num parts',  (num_bytes-8)/float(nparts_var_per_file)
            self._blocks[name] = GadgetBinaryBlock(self, name, elements, self.header.dtype, self._files, nparts_var, file_offsets)
            file_offsets += self.header.dtype_width * nparts_var_per_file * elements + 8 + 16

    
    def __iter__(self):
        if self._blocks is None:
            self._load_blocks()

        return iter(self._blocks.keys())

    def __getitem__(self, name):
        if self._blocks is None:
            self._load_blocks()
        
        return self._blocks[name]

##########################################################################################       
def _get_all_format2_names(filename, swap_endian):
    """ find all the block names in a gadget format-2 file """
    u32 = _rtype(uint32, swap_endian)
    f = open(filename, 'rb')
    blocks = []
    while True:
        eight = f.read(4)
        if len(eight)==0:
            break
        n = np.fromstring(eight, u32,1)[0]
        assert(n==8)
        name = f.read(4)

        blocks.append(name)
        count = np.fromfile(f, u32,1)[0]
        f.read(4)
        count2 = np.fromfile(f, u32,1)[0]
        f.seek(count2+4, 1)
    f.close()
    return blocks




def load_snapshot(filename="", directory="", snapnum=None, additional_blocks=None, label_table=None):
    """
    Returns a Gadget snapshot object. Loads both binary format 1 and format 2 snapshots.    
    Keyword arguments:
        filename:          -- the complete filename or the root of the snapshot filenames
        directory          -- the directory in which the snapshot files or directories can 
                              be found
        snapnum:           -- which snapshot to load
        additional_blocks  -- a list of blocks other than the basic ones that exist in the 
                              snapshot (used for format 1 only)
        label_table:       -- a dictionary of additional labels (e.g. {"ACCE":(3, (1,1,1,1,1,1))})
                              which specifies the number of fields and which particles this label is for.
                              (used for format 2 only)

    """
    
    files = _get_files(directory, filename, snapnum)
    
    format = _determine_binary_format(files[0])
    if format==1:
        return GadgetBinaryFormat1Snapshot(files=files, additional_blocks=additional_blocks)
    elif format==2:
        return GadgetBinaryFormat2Snapshot(files=files, label_table=label_table)
        
def load_snapshot_files(filename="", directory="", snapnum=None, additional_blocks=None):
    files = _get_files(directory, filename, snapnum)
    
    format = _determine_binary_format(files[0])    

    if format==1:
        for file in files:
            yield GadgetBinaryFormat1Snapshot(files=[file], additional_blocks=additional_blocks)
    elif format==2:       
        for file in files:      
            yield GadgetBinaryFormat2Snapshot(files=[file])           

def save_snapshot(filename, header, pos, vel, ids=None, masses=None, extra_data=None):
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
