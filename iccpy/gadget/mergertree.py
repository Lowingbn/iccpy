import h5py
import numpy as np
import iccpy.utils

class MergerTree:
    def __init__(self, filename):
        print "Building merger tree..."

        file = h5py.File(filename, 'r')
        self.nodeID = np.array(file['/haloTrees/nodeIndex'])
        descendantID = np.array(file['/haloTrees/descendantIndex'])
        self.snapshotNum = np.array(file['/haloTrees/snapshotNumber'])
        isMainProgenitor = np.array(file['/haloTrees/isMainProgenitor'])

        num = len(self.nodeID)

        #Construct descendantIndex
        self.descendantIndex = iccpy.utils.match(descendantID, self.nodeID)

        #Construct progentiorIndex
        idx = np.where(isMainProgenitor==1)[0]
        self.mainProgenitorIndex = iccpy.utils.match(self.nodeID, descendantID[idx])
        self.mainProgenitorIndex[np.where(self.mainProgenitorIndex!=-1)] = idx[self.mainProgenitorIndex[np.where(self.mainProgenitorIndex!=-1)]]

    def get_descendants(self, subhalo_id):
        idx = np.where(self.nodeID==subhalo_id)[0][0]

        descendants = []
        while idx!=-1:
            descendants.append(self.nodeID[idx])
            idx = self.descendantIndex[idx]

        return np.array(descendants)

    def get_main_progenitors(self, subhalo_id):
        idx = np.where(self.nodeID==subhalo_id)[0][0]

        progenitors = []
        while idx!=-1:
            progenitors.append(self.nodeID[idx])
            idx = self.mainProgenitorIndex[idx]

        return np.array(progenitors)

    def get_progenitors(self, subhalo_id):
        sub_idx = np.where(self.nodeID==subhalo_id)[0][0]
        
        idx = np.where(self.descendantIndex==sub_idx)
        return self.nodeID[idx]


if __name__ == "__main__":
    #mt = MergerTree("/Users/bjlowin/Work/aquarius/halo_data/Aq-A/5/tree_127.0.hdf5")
    mt = MergerTree("/gpfs/data/jch/Aquarius/Trees/Aq-A/4/trees/treedir_127/tree_127.0.hdf5")

    print mt.get_main_progenitors(32003800000109)
    print mt.get_descendants(32003800000109)
