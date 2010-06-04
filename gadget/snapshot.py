import binary_snapshot_io
import os.path

def get_filename(directory, file, snapnum=0, filenum=0):
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

def load_snapshot_file(directory, filename, snapshot_num, file_num=0):
    filename = get_filename(directory, filename, snapshot_num, file_num)
    return binary_snapshot_io.read_snapshot_file(filename)

def load_snapshot_files(directory, file, snapshot_num):
    filename = get_filename(directory, file, snapshot_num)
    
    header = binary_snapshot_io.read_header(filename)
    num_files = header['num_files']
    
    for i in range(len(num_files)):
        filename = get_filename(directory, file, snapshot_num, i)
        header, res = binary_snapshot_io.read_snapshot_file(filename, False)
        yield header, res
