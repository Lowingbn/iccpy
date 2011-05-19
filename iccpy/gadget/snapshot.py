
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

def convert_to_physical(header, res):
    """ convert data from a Gadget snapshot into physical coordinates """
    h = header['hubble0']
    z = header['redshift']
    
    res['pos'] *= (1+z) / h
    res['vel'] = res['vel'] * np.sqrt(1+z) + iccpy.cosmology.hubble_param(1/(1+z), h) * res['pos']
    res['mass'] /= h
    
    if 'sml' in res: res['sml'] *= (1+z) / h

def convert_to_comoving(header, res):
    """ convert data from a Gadget snapshot into co-moving coordinates """
    h = header['hubble0']
    z = header['redshift']
