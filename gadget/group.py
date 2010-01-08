import binary_group_io

class FoFGroup:
	"""Store information about each FoF group"""
	def __init__(self, snapshot_num):
		self.snapshot_num = snapshot_num
		self.npart = 0
	

class SubfindGroup:
	"""Store information about each subfind group"""
	def __init__(self, snapshot_num):
		self.snapshot_num = snapshot_num
		self.npart = 0
		self.pot_min = None
		self.most_bound_particle = None
		self.mass = 0
		self.half_mass_radius = 0
		self.v_max = 0
		

def load_FoF_groups(directory, snapshot_num, type='bin', loadIDs=False):
	if type=='bin':
		if loadIDs:
			filebase = "%s/groups_%03d/group_ids_%03d" % (directory, snapshot_num, snapshot_num)
			ids = binary_group_io.read_IDs(filebase)
		else:
			ids = None

		filename = "%s/groups_%03d/group_tab_%03d.0" % (directory, snapshot_num, snapshot_num)
		return binary_group_io.read_group_file(filename, ids)
	else:
		pass
			
def load_subfind_groups(directory, snapshot_num, type='bin', loadIDs=False):
	if type=='bin':
		if loadIDs:
			filebase = "%s/groups_%03d/subhalo_ids_%03d" % (directory, snapshot_num, snapshot_num)
			ids = binary_group_io.read_IDs(filebase)
		else:
			ids = None

		filename = "%s/groups_%03d/subhalo_tab_%03d.0" % (directory, snapshot_num, snapshot_num)
		return binary_group_io.read_subfind_file(filename, ids)
	else:
		pass
	
def convert_to_subfind_groups(props, snapshot_num):
	subgroups = []
	
	for i in range(len(props['npart'])):
		subgroup = SubfindGroup(snapshot_num)
		subgroup.npart = props['npart'][i]
		subgroup.pot_min = props['pot_min'][i]
		subgroup.most_bound_particle = props['most_bound_particle'][i]
		subgroup.mass = props['mass'][i]
		subgroup.v_max = props['vel_max'][i]
		subgroup.half_mass_radius = props['half_mass_radius'][i]		
		
		if 'ids' in props:
			subgroup.ids = props['ids'][i]
		
		subgroups.append(subgroup)
	
	return subgroups