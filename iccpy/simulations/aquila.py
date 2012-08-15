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
    print group

    header = Struct()
    return header
    
if __name__=="__main__":
    header = read_header("test")
    print header.a
