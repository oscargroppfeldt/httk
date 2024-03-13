




#from httk.atomistic import Structure
#structure = Structure.io.load('/home/abhijith/work/httk/Tutorial/Step2/POSCAR')
#print(structure.pc)

#from httk.iface.vasp_if import structure_to_poscar
#%%
from httk.atomistic import Structure
structure = Structure.io.load('/home/abhijith/work/httk/Tutorial/Step2/POSCAR')
transformation=[[3,0,0],[0,3,0],[0,0,3]]
st=structure.transform(transformation,pass_sym=True)
#structure_to_poscar('POSCAR_PC', st,fix_negative_determinant=True, comment='vanilla')
#%%

'''
import os, errno

import httk
import httk.db
from httk.atomistic import Structure

# Remove example database
try:
    os.remove('st___.sqlite')
except OSError as e:
    if e.errno != errno.ENOENT:
        raise

backend = httk.db.backend.Sqlite('st___.sqlite')
store = httk.db.store.SqlStore(backend)

struct = Structure.io.load('/home/abhijith/work/httk/Tutorial/Step2/POSCAR')

transformation=[[3,0,0],[0,3,0],[0,0,3]]
st=struct.transform(transformation,pass_sym=True)
store.save(st)

search = store.searcher()
search_struct = search.variable(Structure)
#search.add(search_struct.uc_nbr_atoms < 40)

search.output(search_struct, 'structure ')

for match in search:
    structure = match[0][0]
    print("Found:", struct.formula)

from httk.iface.vasp_if import structure_to_poscar



from httk.atomistic import Structure
structure = Structure.io.load('/home/abhijith/work/httk/Tutorial/Step2/POSCAR')
#print(structure.pc)
print(structure.uc.uc_sites.counts)
transformation=[[3,0,0],[0,3,0],[0,0,3]]
st_2=structure.transform(transformation,pass_sym=False)
#print(st_2.uc.uc_sites.counts)
#st_2.io.save('st_after_mod_uc.vasp')
#print(st_2.rc.rc_sites.counts)

structure_to_poscar('st_POSCAR_ab', st_2,fix_negative_determinant=True, comment='abhijith')




from httk.atomistic import Structure
structure = Structure.io.load('/home/abhijith/work/httk/Tutorial/Step2/POSCAR')
transformation=[[3,0,0],[0,3,0],[0,0,3]]
st_2=structure.transform(transformation,pass_sym=False)
print('first print of uc_basis: \n', st_2.uc_basis)
print(st_2.uc.uc_sites.counts)
print('representations: \n',st_2._other_reps)
print(st_2.pc.uc_basis)
print(st_2._other_reps)
print('second print of uc_basis: \n',st_2.uc_basis)
print(st_2.uc.uc_sites.counts)



from httk.atomistic import Structure
structure = Structure.io.load('/home/abhijith/work/httk/Tutorial/Step2/POSCAR')
print(structure.uc)
print(structure.pc)
print(structure.uc)


from httk.atomistic import Structure
structure = Structure.io.load('/home/abhijith/work/httk/Tutorial/Step2/POSCAR')
print(structure.uc)
print(structure.rc_cell.get_axes_standard_order_transform())


from httk.atomistic import Structure
from httk.atomistic.structureutils import get_primitive_basis_transform
structure = Structure.io.load('/home/abhijith/work/httk/Tutorial/Step2/POSCAR')
print(structure.uc)
print(get_primitive_basis_transform(structure.rc_sites.hall_symbol))
print(structure.rc_cell.get_axes_standard_order_transform())


# construct defect cell
final_assignments_class = Assignments.create(assignments=temp_assignments)
rc_sites = RepresentativeSites.create(reduced_coordgroups=final_coordgroups, spacegroup='P 1')
#rc_sites = RepresentativeSites(reduced_coordgroups=final_coordgroups, hall_symbol='P 1')
final_struct = Structure(final_assignments_class, rc_ites=rc_sites, rc_cell=struct.rc_cell)




ase_stru=see(structure)

ase_spr=make_supercell(ase_stru,transformation)

st=structure.transform(transformation)

ase_stru_supr=see(st)


import cProfile
cProfile.run('st=structure.transform(transformation)',sort='time',filename='/home/abhijith/work/log_folder/httk_supercell')

import pstats
p=pstats.Stats('/home/abhijith/work/log_folder/httk_supercell')
p.sort_stats('cumtime').print_stats(20)



import ase.build.tools as tool

ase_stru_supr_sort=tool.sort(ase_stru_supr)
ase_spr_sort=tool.sort(ase_spr)


from ase.utils.structure_comparator import SymmetryEquivalenceCheck
SymmetryEquivalenceCheck().compare(ase_stru_supr,ase_spr_sort)


import ase.spacegroup.spacegroup as sp
spgp=sp.get_spacegroup(ase_stru)


structure = Structure.io.load("/home/abhijith/work/joel_work/ADAQ/2_Unitcell_Workflow/Runs/ht.task.tetralith--default.diamond_pbe.cleanup.0.unclaimed.3.finished/ht.run.2021-06-29_16.12.18/CONTCAR.relaxfinal.bz2")


from httk.atomistic.spacegrouputils import get_symops_strs, spacegroup_get_number_and_setting

spacegroup_get_number_and_setting('F 4d 2 3 -1d')




struct=structure


import httk.atomistic.data
from httk.core.httkobject import HttkPlugin, HttkPluginWrapper
from httk.atomistic import Cell
from httk.atomistic.spacegrouputils import get_symops_strs, spacegroup_get_number_and_setting

from httk import config
from httk.atomistic import Structure, UnitcellSites
import httk.iface
from httk.iface.ase_if import *
from httk.external.subimport import submodule_import_external
import ase
import ase.io
import ase.geometry
from ase.spacegroup import crystal
from ase.atoms import Atoms

try:
	from ase import version

	ase_version = version.version
except ImportError:
	from ase import __version__ as aseversion

ase_major_version = aseversion.split('.')[0]
ase_minor_version = aseversion.split('.')[1]

if struct.has_uc_repr:
	symbollist, scaled_positions = httk.iface.ase_if.uc_structure_to_symbols_and_scaled_positions(struct)
	scaled_positions = scaled_positions.to_floats()
	cell = struct.uc_basis.to_floats()
	symbols = []
	for s in symbollist:
		if is_sequence(s):
			if len(s) == 1:
				symbols += [s[0]]
			else:
				symbols += [str(s)]
		else:
			symbols += [s]

	atoms = Atoms(symbols=symbols,
	              cell=cell,
	              scaled_positions=scaled_positions,
	              pbc=True)

elif struct.has_rc_repr:
	symbollist, scaled_positions = httk.iface.ase_if.rc_structure_to_symbols_and_scaled_positions(struct)
	symbols = []
	for s in symbollist:
		if is_sequence(s):
			if len(s) == 1:
				symbols += [s[0]]
			else:
				symbols += [str(s)]
		else:
			symbols += [s]

	hall = struct.rc_sites.hall_symbol
	symops = get_symops_strs(hall)
	rot, trans = ase.spacegroup.spacegroup.parse_sitesym(symops)
	spgnbr, setting = spacegroup_get_number_and_setting(hall)
	spg = ase.spacegroup.spacegroup.spacegroup_from_data(no=spgnbr,
	                                                     symbol=None,
	                                                     setting=int(setting),
	                                                     centrosymmetric=None,
	                                                     scaled_primitive_cell=None,
	                                                     reciprocal_cell=None,
	                                                     subtrans=None,
	                                                     sitesym=symops,
	                                                     rotations=rot,
	                                                     translations=trans,
	                                                     datafile=None)  # ase not reading hall number correctly, changed to spgnbr and setting,symops

	atoms = crystal(symbols=symbols, basis=scaled_positions.to_floats(), spacegroup=spg,
	                cellpar=[float(struct.rc_a), float(struct.rc_b), float(struct.rc_c),
	                         float(struct.rc_alpha), float(struct.rc_beta), float(struct.rc_gamma)])

'''