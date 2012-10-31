import iccpy.gadget

sim_label = { 'aqa' : 'A', 'aqb' : 'B', 'aqc':'C', 'aqd':'D', 'aqe':'E' }
last_snapnum = { 'aqa2' : 1023, 'aqa3' : 511, 'aqa4' : 1023, 'aqb2' : 127, 'aqc2' : 127, 'aqd2' : 127, 'aqe2' : 127 }

r_200 = { 'aqa1' : 245.67, 'aqa2' : 245.88, 'aqa3' : 245.64, 'aqa4' : 245.70, 'aqa5' : 246.37, 'aqb2' : 187.70, 'aqb4' : 188.85,
          'aqc2' : 242.82, 'aqc4' : 243.68, 'aqd2' : 242.85, 'aqd4' : 243.60, 'aqe2' : 212.28, 'aqe4' : 213.63, 'aqf2' : 209.21,
          'aqf4' : 207.15 }
M_200 = { 'aqa1' : 183.9,  'aqa2' : 184.2,  'aqa3' : 183.6,  'aqa4' : 183.8,  'aqa5' : 185.3,  'aqb2' : 81.94,  'aqb4' : 83.45,
          'aqc2' : 177.4,  'aqc4' : 179.3,  'aqd2' : 177.4,  'aqd4' : 179.1,  'aqe2' : 118.5,  'aqe4' : 120.8,  'aqf2' : 113.5,
          'aqf4' : 110.1 }

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
    print load_last_snapshot("aqa4")
    print get_halo_centre("aqa4")
