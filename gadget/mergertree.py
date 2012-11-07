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
        self.mass = np.array(file['/haloTrees/nodeMass'])

        num = len(self.nodeID)

        #Map snapshots to expansion factors
        all_snapshot_nums = np.array(file['/outputTimes/snapshotNumber'])
        exp_factors = np.array(file['/outputTimes/expansion'])

        self.exp_factors = np.zeros(np.max(all_snapshot_nums)+1)
        self.exp_factors[all_snapshot_nums] = exp_factors

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

    def get_max_mass(self, subhalo_id):
        idx = np.where(self.nodeID==subhalo_id)[0][0]
        prog_idxs = []
        while idx!=-1:
            prog_idxs.append(idx)
            idx = self.mainProgenitorIndex[idx]

        max_mass_idx = np.argmax(self.mass[prog_idxs])
        return self.mass[prog_idxs[max_mass_idx]], self.exp_factors[self.snapshotNum[prog_idxs[max_mass_idx]]+1]


if __name__ == "__main__":
    #mt = MergerTree("/Users/bjlowin/Work/aquarius/halo_data/Aq-A/5/tree_127.0.hdf5")
    mt = MergerTree("/gpfs/data/jch/Aquarius/Trees/Aq-A/4/trees/treedir_127/tree_127.0.hdf5")
    #mt = MergerTree("/gpfs/data/jch/Aquarius/Trees/Aq-A/2/trees/treedir_127/tree_127.0.hdf5")

    print mt.get_main_progenitors(127000000000000)
    print mt.get_descendants(127000000000000)
    print mt.get_max_mass(127000000000000)
