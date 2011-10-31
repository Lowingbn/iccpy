
import binary_snapshot_io
import os.path
import numpy as np
import iccpy.cosmology

def _get_filename(directory, file, snapnum=0, filenum=0):
    if filenum==0 and snapnum==0:
        filename = directory + "/" + file
        if os.path.isfile(filename): return filename 
    
    if filenum==0:
        filename = "%s/%s_%03d" % (directory, file, snapnum)
        if os.path.isfile(filename): return filename
        
        filename = "%s/snapdir_%03d/%s_%03d" % (directory, snapnum, file, snapnum)
        if os.path.isfile(filename): return filename
        
    filename = "%s/%s_%03d.%d" % (directory, file, snapnum, filenum)
    if os.path.isfile(filename): return filename
    
    filename = "%s/snapdir_%03d/%s_%03d.%d" % (directory, snapnum, file, snapnum, filenum)
    if os.path.isfile(filename): return filename
    
    print "Unable to find Gadget snapshot file %s %s %d %d.\n" % (directory, file, snapnum, filenum)
    raise IOError
  
def load_header(directory, filename, snapshot_num, file_num=0):
    return load_snapshot_header(directory, filename, snapshot_num, file_num)
    
def load_snapshot_header(directory, filename, snapshot_num, file_num=0):
    """ load a single Gadget snapshot file header """
    filename = _get_filename(directory, filename, snapshot_num, file_num)
    return binary_snapshot_io.read_snapshot_header(filename)

def load_file(directory, filename, snapshot_num, file_num=0):
    return load_snapshot_file(directory, filename, snapshot_num, file_num)

def load_snapshot_file(directory, filename, snapshot_num, file_num=0):
    """ load a single Gadget snapshot file """
    filename = _get_filename(directory, filename, snapshot_num, file_num)
    return binary_snapshot_io.read_snapshot_file(filename)

def load(directory, file, snapshot_num):
    return load_snapshot_files(directory, file, snapshot_num)

def load_snapshot_files(directory, file, snapshot_num):
    """Loads a complete binary GADGET snapshot
    
    Arguments:
    directory -- the location of the snapshot
    file -- the file name
    snapshot_num -- the number of the snapshot 
    """
    filename = _get_filename(directory, file, snapshot_num)
    
    header = binary_snapshot_io.read_snapshot_header(filename)
    num_files = header['num_files']
    
    for i in range(num_files[0]):
        filename = _get_filename(directory, file, snapshot_num, i)
        header, res = binary_snapshot_io.read_snapshot_file(filename, False)
        yield header, res
        
def load_snapshot(directory, file, snapshot_num):
    """Combines a set of snapshot files to load a complete binary GADGET snapshot
    
    Arguments:
    directory -- the location of the snapshot
    file -- the file name
    snapshot_num -- the number of the snapshot 
    """
    filename = _get_filename(directory, file, snapshot_num)
    header, res = binary_snapshot_io.read_snapshot_file(filename, False)
    num_files = header['num_files'][0]
    
    pos = [ np.empty([header['npartTotal'][i], 3]) for i in range(6)]
    vel = [ np.empty([header['npartTotal'][i], 3]) for i in range(6)]
    ids = [ np.empty(header['npartTotal'][i], dtype=res['id'].dtype) for i in range(6)]
    mass = [ 0 if header['mass'][i]==0 else np.empty(header['npartTotal'][i]) for i in range(6)]
    count = np.zeros(6)

    for i in range(num_files):
        filename = _get_filename(directory, file, snapshot_num, i)
        h, res = binary_snapshot_io.read_snapshot_file(filename, False)
        
        idxs = np.insert(np.cumsum(h['npart']), 0, 0)        
        for j in range(6):
            pos[j][count[j]:count[j]+h['npart'][j]] = res['pos'][idxs[j]:idxs[j+1]]
            vel[j][count[j]:count[j]+h['npart'][j]] = res['vel'][idxs[j]:idxs[j+1]]
            ids[j][count[j]:count[j]+h['npart'][j]] = res['id'][idxs[j]:idxs[j+1]]
            if header['mass'][i]==0:
                mass[j][count[j]:count[j]+h['npart'][j]] = res['mass'][idxs[j]:idxs[j+1]]
            
            count[j] += h['npart'][j]
        
    header['num_files'][0] = 1
    header['npart'] = header['npartTotal']
    
    res['pos'] = pos
    res['vel'] = vel
    res['mass'] = mass
    res['id'] = ids
    
    return header, mass


def convert_to_physical(header, res):
    """ convert data from a Gadget snapshot into physical coordinates """
    z = header['redshift']
    a = 1/(1+z)
    
    res['pos'] *= a
    res['vel'] = res['vel'] * np.sqrt(a) + iccpy.cosmology.hubble_param(a) * res['pos']
    
    if 'sml' in res: res['sml'] *= (1+z)

def convert_to_comoving(header, res):
    """ convert data from a Gadget snapshot into co-moving coordinates """
    h = header['hubble0']
    z = header['redshift']
