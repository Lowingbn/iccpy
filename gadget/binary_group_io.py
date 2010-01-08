import os.path
import numpy as np
from numpy import uint32, uint64, float64, float32

def read_group_file(filename, ids=None):
	""" Reads a binary FoF group file """
	f = open(filename, mode='rb')

	num_groups = np.fromfile(f, uint32, 1)[0]
	num_groups_tot = np.fromfile(f, uint32, 1)[0]
	num_ids = np.fromfile(f, uint32, 1)[0]
	num_ids_tot = np.fromfile(f, uint64, 1)[0]
	num_files = np.fromfile(f, uint32, 1)[0]

	res = {}
	
	res['npart'] = np.fromfile(f, uint32, num_groups)
	offsets = np.fromfile(f, uint32, num_groups)
	res['mass'] = np.fromfile(f, float32, num_groups)
	res['com'] = np.fromfile(f, float32, 3*num_groups).reshape((num_groups, 3))
	res['vel'] = np.fromfile(f, float32, 3*num_groups).reshape((num_groups, 3))

	res['typeLen'] = np.fromfile(f, uint32, 6*num_groups).reshape((num_groups, 6))
	res['typeMass'] = np.fromfile(f, float32, 6*num_groups).reshape((num_groups, 6))
  
	if ids is not None:
		res['id'] = []
		for offset, len in zip(offsets, res['npart']):
			res['id'].append(ids[offset:offset+len])

	f.close()
	
	return num_groups, res
	
def read_subfind_file(filename, ids=None):
	""" Reads a binary Subfind group file """
	f = open(filename, mode='rb')
	
	num_groups = np.fromfile(f, uint32, 1)[0]
	num_groups_tot = np.fromfile(f, uint32, 1)[0]
	num_ids = np.fromfile(f, uint32, 1)[0]
	num_ids_tot = np.fromfile(f, uint64, 1)[0]
	num_files = np.fromfile(f, uint32, 1)[0]
	num_subgroups = np.fromfile(f, uint32, 1)[0]
	num_subgroups_tot = np.fromfile(f, uint32, 1)[0]	
	
	groups = {}
	
	groups['npart'] = np.fromfile(f, uint32, num_groups)
	groups['offsets'] = np.fromfile(f, uint32, num_groups)
	groups['mass'] = np.fromfile(f, float32, num_groups)
	groups['pot_min'] = np.fromfile(f, float32, 3*num_groups).reshape((num_groups, 3))
	
	groups['mass_mean_200'] = np.fromfile(f, float32, num_groups)
	groups['radius_mean_200'] = np.fromfile(f, float32, num_groups)
	groups['mass_crit_200'] = np.fromfile(f, float32, num_groups)
	groups['radius_crit_200'] = np.fromfile(f, float32, num_groups)
	groups['mass_top_hat_200'] = np.fromfile(f, float32, num_groups)
	groups['radius_top_hat_200'] = np.fromfile(f, float32, num_groups)
	
	groups['npart_contaimination'] = np.fromfile(f, uint32, num_groups)
	groups['mass_contaimination'] = np.fromfile(f, uint32, num_groups)
	
	groups['num_subgroups'] = np.fromfile(f, uint32, num_groups)  
	groups['subgroup_offset'] = np.fromfile(f, uint32, num_groups)
	
	if ids is not None:
		groups['ids'] = []
		for offset, len in zip(groups['offsets'], groups['npart']):
			groups['ids'].append(ids[offset:offset+len])
	
	subgroups = {}
	
	group_idx = [ [i]*x for i,x in enumerate(groups['num_subgroups']) ]
	subgroups['group_idx'] = np.concatenate(group_idx)
	
	subgroups['npart'] = np.fromfile(f, uint32, num_subgroups)
	subgroups['offsets'] = np.fromfile(f, uint32, num_subgroups)
	subgroups['parent'] = np.fromfile(f, uint32, num_subgroups)
	subgroups['mass'] = np.fromfile(f, float32, num_subgroups)
	subgroups['pot_min'] = np.fromfile(f, float32, 3*num_subgroups).reshape((num_subgroups, 3))
	subgroups['vel'] = np.fromfile(f, float32, 3*num_subgroups).reshape((num_subgroups, 3))
	subgroups['com'] = np.fromfile(f, float32, 3*num_subgroups).reshape((num_subgroups, 3))
	subgroups['spin'] = np.fromfile(f, float32, 3*num_subgroups).reshape((num_subgroups, 3))
	
	subgroups['vel_disp'] = np.fromfile(f, float32, num_subgroups)
	subgroups['vel_max'] = np.fromfile(f, float32, num_subgroups)
	subgroups['radius_vel_max'] = np.fromfile(f, float32, num_subgroups)
	subgroups['half_mass_radius'] = np.fromfile(f, float32, num_subgroups)
	
	#From length of file we can work out size of id
	filelength = os.path.getsize(filename)	
	if filelength==32 + 16*4*num_groups + 23*4*num_subgroups: idType = uint64
	elif filelength==32 + 16*4*num_groups + 22*4*num_subgroups: idType = uint32
	else: raise Exception('Unable to determine size of ID type from file length')
	
	subgroups['most_bound_particle'] = np.fromfile(f, idType, num_subgroups)
	subgroups['subgroups'] = np.fromfile(f, uint32, num_subgroups)

	if ids is not None:
		subgroups['ids'] = []
		for offset, len in zip(groups['offsets'], subgroups['npart']):
			subgroups['ids'].append(ids[offset:offset+len])
				
	return num_groups, groups, num_subgroups, subgroups
	
def read_IDs(filebase):
	""" Reads a binary ID group file """
	filename = "%s.0" % (filebase)
	f = open(filename, mode='rb')
	f.seek(8, 1)
	num_ids = np.fromfile(f, uint32, 1)[0]
	num_ids_total = np.fromfile(f, uint64, 1)[0]
	num_files = np.fromfile(f, uint32, 1)[0]
	f.close()
	
	filelength = os.path.getsize(filename)
	
	if filelength==28 + 8*num_ids: idType = uint64
	elif filelength==28 + 4*num_ids: idType = uint32
	else: raise Exception('Unable to determine size of ID type from file length')
	
	idList = []
	
	for i in range(num_files):
		filename = "%s.%d" % (filebase, i)
		f = open(filename, mode='rb')
		
		f.seek(8, 1)
		num_ids = np.fromfile(f, uint32, 1)[0]
		f.seek(16, 1)
		ids = np.fromfile(f, idType, num_ids)
		f.close()
	
		idList.append(ids)
	
	ids = np.concatenate(idList)

	return ids