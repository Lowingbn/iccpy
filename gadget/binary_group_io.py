import os.path
import numpy as np

def read_group_file(filename, ids=None, swap_endian=False):
	""" Reads a binary FoF group file """
	f = open(filename, mode='rb')
	
	if swap_endian:
		u32 = 'u4'
		f32 = 'f4'
	else:
		u32 = '>u4'
		f32 = '>f4'

	num_groups = np.fromfile(f, u32, 1)[0]
	num_groups_tot = np.fromfile(f, u32, 1)[0]
	num_ids = np.fromfile(f, u32, 1)[0]
	num_ids_tot = np.fromfile(f, f32, 1)[0]
	num_files = np.fromfile(f, u32, 1)[0]
	
	props = {}
	props['num_groups'] = num_groups
	props['num_groups_tot'] = num_groups_tot
	props['num_ids'] = num_ids
	props['num_ids_tot'] = num_ids_tot
	props['num_files'] = num_files
	
	res = {}
	
	res['npart'] = np.fromfile(f, u32, num_groups)
	offsets = np.fromfile(f, u32, num_groups)
	res['mass'] = np.fromfile(f, f32, num_groups)
	res['com'] = np.fromfile(f, f32, 3 * num_groups).reshape((num_groups, 3))
	res['vel'] = np.fromfile(f, f32, 3 * num_groups).reshape((num_groups, 3))

	res['typeLen'] = np.fromfile(f, u32, 6 * num_groups).reshape((num_groups, 6))
	res['typeMass'] = np.fromfile(f, f32, 6 * num_groups).reshape((num_groups, 6))
	
	if ids is not None:
		res['id'] = []
		for offset, len in zip(offsets, res['npart']):
			res['id'].append(ids[offset:offset + len])

	f.close()
	
	return props, res
	
def read_subfind_file(filename, ids=None, swap_endian=None):
	""" Reads a binary Subfind group file """
	f = open(filename, mode='rb')
	
	if swap_endian:
		u32 = 'u4'
		f32 = 'f4'
		u64 = 'u8'
	else:
		u32 = '>u4'
		f32 = '>f4'
		u64 = '>u8'
	
	num_groups = np.fromfile(f, u32, 1)[0]
	num_groups_tot = np.fromfile(f, u32, 1)[0]
	num_ids = np.fromfile(f, u32, 1)[0]
	num_ids_tot = np.fromfile(f, u64, 1)[0]
	num_files = np.fromfile(f, u32, 1)[0]
	num_subgroups = np.fromfile(f, u32, 1)[0]
	num_subgroups_tot = np.fromfile(f, u32, 1)[0]
	
	props = {}
	props['num_groups'] = num_groups
	props['num_groups_tot'] = num_groups_tot
	props['num_ids'] = num_ids
	props['num_ids_tot'] = num_ids_tot
	props['num_files'] = num_files
	props['num_subgroups'] = num_subgroups
	props['num_subgroups_tot'] = num_subgroups_tot	
	
	groups = {}
	
	groups['npart'] = np.fromfile(f, u32, num_groups)
	groups['offsets'] = np.fromfile(f, u32, num_groups)
	groups['mass'] = np.fromfile(f, f32, num_groups)
	groups['pot_min'] = np.fromfile(f, f32, 3 * num_groups).reshape((num_groups, 3))
	
	groups['mass_mean_200'] = np.fromfile(f, f32, num_groups)
	groups['radius_mean_200'] = np.fromfile(f, f32, num_groups)
	groups['mass_crit_200'] = np.fromfile(f, f32, num_groups)
	groups['radius_crit_200'] = np.fromfile(f, f32, num_groups)
	groups['mass_top_hat_200'] = np.fromfile(f, f32, num_groups)
	groups['radius_top_hat_200'] = np.fromfile(f, f32, num_groups)
	
	groups['npart_contaimination'] = np.fromfile(f, u32, num_groups)
	groups['mass_contaimination'] = np.fromfile(f, u32, num_groups)
	
	groups['num_subgroups'] = np.fromfile(f, u32, num_groups)  
	groups['subgroup_offset'] = np.fromfile(f, u32, num_groups)
	
	if ids is not None:
		groups['ids'] = []
		for offset, length in zip(groups['offsets'], groups['npart']):
			groups['ids'].append(ids[offset:offset + length])
	
	subgroups = {}
	
	subgroups['npart'] = np.fromfile(f, u32, num_subgroups)
	subgroups['offsets'] = np.fromfile(f, u32, num_subgroups)
	subgroups['parent'] = np.fromfile(f, u32, num_subgroups)
	subgroups['mass'] = np.fromfile(f, f32, num_subgroups)
	subgroups['pot_min'] = np.fromfile(f, f32, 3 * num_subgroups).reshape((num_subgroups, 3))
	subgroups['vel'] = np.fromfile(f, f32, 3 * num_subgroups).reshape((num_subgroups, 3))
	subgroups['com'] = np.fromfile(f, f32, 3 * num_subgroups).reshape((num_subgroups, 3))
	subgroups['spin'] = np.fromfile(f, f32, 3 * num_subgroups).reshape((num_subgroups, 3))
	
	subgroups['vel_disp'] = np.fromfile(f, f32, num_subgroups)
	subgroups['vel_max'] = np.fromfile(f, f32, num_subgroups)
	subgroups['radius_vel_max'] = np.fromfile(f, f32, num_subgroups)
	subgroups['half_mass_radius'] = np.fromfile(f, f32, num_subgroups)
	
	#From length of file we can work out size of id
	filelength = os.path.getsize(filename)	
	if filelength == 32 + 16 * 4 * num_groups + 23 * 4 * num_subgroups: idType = u64
	elif filelength == 32 + 16 * 4 * num_groups + 22 * 4 * num_subgroups: idType = u32
	else: raise Exception('Unable to determine size of ID type from file length')
	
	subgroups['most_bound_particle'] = np.fromfile(f, idType, num_subgroups)
	subgroups['subgroups'] = np.fromfile(f, u32, num_subgroups)

	if ids is not None:
		subgroups['ids'] = []
		for offset, length in zip(subgroups['offsets'], subgroups['npart']):
			subgroups['ids'].append(ids[offset:offset + length])
							
	return props, groups, subgroups
	
def read_IDs(filebase, swap_endian=False):
	""" Reads a binary ID group file """
	if swap_endian:
		u32 = 'u4'
		u64 = 'u8'
	else:
		u32 = '>u4'
		u64 = '>u8'
		
	filename = "%s.0" % (filebase)
	f = open(filename, mode='rb')
	f.seek(8, 1)
	num_ids = np.fromfile(f, u32, 1)[0]
	num_ids_total = np.fromfile(f, u64, 1)[0]
	num_files = np.fromfile(f, u32, 1)[0]
	f.close()
	
	filelength = os.path.getsize(filename)
	
	if filelength == 28 + 8 * num_ids: idType = u64
	elif filelength == 28 + 4 * num_ids: idType = u32
	else: raise Exception('Unable to determine size of ID type from file length')
	
	idList = []
	
	for i in range(num_files):
		filename = "%s.%d" % (filebase, i)
		f = open(filename, mode='rb')
		
		f.seek(8, 1)
		num_ids = np.fromfile(f, u32, 1)[0]
		f.seek(16, 1)
		ids = np.fromfile(f, idType, num_ids)
		f.close()
	
		idList.append(ids)
	
	ids = np.concatenate(idList)

	return ids
