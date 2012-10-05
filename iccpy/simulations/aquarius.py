import iccpy.gadget

sim_label = { 'aqa' : 'A', 'aqb' : 'B', 'aqc':'C', 'aqd':'D', 'aqe':'E' }
last_snapnum = { 'aqa2' : 1023, 'aqa3' : 511, 'aqa4' : 1023 }

def get_dir(sim_name):
    return "/gpfs/data/aquarius/halo_data/Aq-%s/%c/" % (sim_label[sim_name[0:3]], sim_name[3])

def load_last_snapshot(sim_name):
    return iccpy.gadget.Snapshot(directory=get_dir(sim_name), snapnum=last_snapnum[sim_name])

def get_subhaloes(sim_name, snapnum=None):
    if snapnum==None:
        snapnum=last_snapnum[sim_name]

    catalogue = iccpy.gadget.SubfindCatalogue(get_dir(sim_name), snapnum)
    return catalogue.subhalo

def get_halo_centre(sim_name):
    return get_subhaloes(sim_name)[0].pot_min

if __name__=="__main__":
    print load_last_snapshot("aqa2")
    print get_halo_centre("aqa2")
