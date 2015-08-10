
from prody.proteins import parsePDB, writePDB
from prody.trajectory import writePSF
import os.path as op
from collections import namedtuple
import matrix
import math
import numpy as np
from parsebonds import parsePDBConnect, writePDBConnect
"""
# strands axii assumed in the X direction

1. load bases pdbs
2. Create 5' end, 3' end, and internal versions of each pdb
2. construct a list of prody.AtomGroups for the strand
3. Apply the axial translation
4. Apply twist per base (rotate X)
5. Concatenate AtomGroups 
6. Install Phosphate bonds and Hydrogen bonds since we know the atom serial 
    number offsets in the atom group
7. Flip for reverse strand (rotate Y 180) and align translate to the forware strand
"""
dna_pdb_files = [   'A_single.pdb',
                    'C_single.pdb', 
                    'G_single.pdb', 
                    'T_single.pdb'
                ]
LOCAL_DIR = op.dirname(op.abspath(__file__))
DATA_DIR = op.join(LOCAL_DIR, 'data')

"""
the index bounds not including the 5' OH and the 3' H group
in the reference PDB files, where indexing begins at 1
so (3, y) skips an H and an O in atoms 1 and 2 and there's
a capping HCAP Hydrogen in the y + 1 slot
"""
ID_OFFSETS = [ (3, 34),
                (3, 32),
                (3, 35),
                (3, 34)
            ]

# convert to zero-based indexing
IDX_OFFSETS = [(x - 1, y - 1) for x, y in ID_OFFSETS]

# offset of the 5' O atom the 3' Phosphate bonds to 
O3_ID_OFFSET = 11
O3_IDX_OFFSET = O3_ID_OFFSET - 1

O3_ID_OFFSET_INT = 9
O3_IDX_OFFSET_INT = O3_ID_OFFSET_INT - 1

AGBase = namedtuple('AGBase', ['idx', 'ag'])

DELTA_X = 3.38      # distance in Angstroms between bases
RADIUS_Y = 10.175   # Approximate radius of each base in Angstroms, unused

BASE_LUT = np.zeros(128, dtype=int)
BASE_LUT[b'A'[0]] = 0
BASE_LUT[b'C'[0]] = 1
BASE_LUT[b'G'[0]] = 2
BASE_LUT[b'T'[0]] = 3
BASE_LUT[b'a'[0]] = 0
BASE_LUT[b'c'[0]] = 1
BASE_LUT[b'g'[0]] = 2
BASE_LUT[b't'[0]] = 3

def loadBasePDBs():
    base_atom_groups = []
    base_bond_groups = []
    for fn in dna_pdb_files:
        the_file = op.join(DATA_DIR, fn)
        ag = parsePDB(the_file)
        # print(ag.getDataLabels())
        base_atom_groups.append(ag)

        bg = parsePDBConnect(the_file)
        base_bond_groups.append(bg)
    return base_atom_groups, base_bond_groups

base_atom_groups, base_bond_groups = loadBasePDBs()

def create5PrimeAGs(base_atom_groups, base_bond_groups):
    """ Drop the trailing HCAP Hydrogen
    """
    five_p_atom_groups = []
    bgs = []
    for i, ag in enumerate(base_atom_groups):
        ag_copy = ag.copy()
        ag_copy._n_atoms = ag_copy._n_atoms - 1

        new_data = {}
        end_idx = ID_OFFSETS[i][1]

        for key, val in ag_copy._data.items():
            # print(key, len(val[:end_idx]), ag_copy._n_atoms)
            new_data[key] = val[:end_idx]
        
        ag_copy._data = new_data
        new_coords = ag_copy._getCoords()[:end_idx]
        ag_copy._coords = None
        ag_copy.setCoords(new_coords)
        # print("coord check %d (%d, %d) %d" %(ag_copy._n_atoms, len(ag._coords[0]), len(ag_copy._coords[0]), end_idx))
        five_p_atom_groups.append(ag_copy)
        bgs.append([x.copy() for x in base_bond_groups[i][:end_idx]])
    return five_p_atom_groups, bgs
# end def
five_p_atom_groups, five_p_bgs = create5PrimeAGs(base_atom_groups, base_bond_groups)

def create3PrimeAGs(base_atom_groups, base_bond_groups):
    """ Drop the lead OH atoms
    """
    three_p_atom_groups = []
    bgs = []
    for i, ag in enumerate(base_atom_groups):
        ag_copy = ag.copy()
        ag_copy._n_atoms = ag_copy._n_atoms - 2

        new_data = {}
        start_idx = IDX_OFFSETS[i][0]
        
        for key, val in ag_copy._data.items():
            if key == 'serial': # offset the serial
                new_data[key] = val[start_idx:] - 2
            else:
                new_data[key] = val[start_idx:]
        ag_copy._data = new_data

        new_coords = ag_copy._getCoords()[start_idx:]
        ag_copy._coords = None
        ag_copy.setCoords(new_coords)
        
        three_p_atom_groups.append(ag_copy)
        bgs.append([x - 2 for x in base_bond_groups[i][start_idx:]])
    return three_p_atom_groups, bgs
# end def
three_p_atom_groups, three_p_bgs = create3PrimeAGs(base_atom_groups, base_bond_groups)

def createINTAGs(base_atom_groups, base_bond_groups):
    """ Drop the lead OH atoms and trailing HCAP Hydrogen
    For internal bases
    """
    int_atom_groups = []
    bgs = []
    for i, ag in enumerate(base_atom_groups):
        ag_copy = ag.copy()
        ag_copy._n_atoms = ag_copy._n_atoms - 3

        new_data = {}
        start_idx = IDX_OFFSETS[i][0]
        end_idx = ID_OFFSETS[i][1]
        for key, val in ag_copy._data.items():
            if key == 'serial': # offset the serial
                new_data[key] = val[start_idx:end_idx] - 2
            else:
                new_data[key] = val[start_idx:end_idx]
        ag_copy._data = new_data

        new_coords = ag_copy._getCoords()[start_idx:end_idx]
        ag_copy._coords = None
        ag_copy.setCoords(new_coords)

        int_atom_groups.append(ag_copy)
        bgs.append([x - 2 for x in base_bond_groups[i][start_idx:end_idx]])
    return int_atom_groups, bgs
# end def
int_atom_groups, int_bgs = createINTAGs(base_atom_groups, base_bond_groups)

def createStrand(seq, 
                origin, 
                name="strand", 
                is_fwd=False, 
                bases_per_turn=10.5,
                theta_offset=0.0,
                create_psf=False):
    """
    """
    # global base_atom_groups
    twist_per_segment = 2.*math.pi/bases_per_turn

    if isinstance(seq, str):
        seq = seq.encode('utf-8')
    ag_list = []
    bg_list = []
    lim = len(seq) - 1
    for i, base in enumerate(seq):
        idx = BASE_LUT[base]
        if i == 0:
            ag = five_p_atom_groups[idx].copy()
            bg = five_p_bgs[idx].copy()
        elif i < lim:
            ag = int_atom_groups[idx].copy()
            bg = int_bgs[idx].copy()

        else:
            ag = three_p_atom_groups[idx].copy()
            bg = three_p_bgs[idx].copy()
        m = matrix.makeTranslation(i*DELTA_X, 0, 0)
        new_coords = matrix.applyTransform(ag._getCoords(), m)
        m = matrix.makeRotationX(i*twist_per_segment)
        new_coords = matrix.applyTransform(new_coords, m)

        ag.setCoords(new_coords)
        ag_list.append(AGBase(idx, ag))
        bg_list.append(bg)

    for x in ag_list:
        print(x)

    # prody.AtomGroup doesn't implement __radd__ so we can't do sum
    ag_out = ag_list[0].ag
    bg_out = bg_list[0]
    offset = ag_out.numAtoms()
    offset__idx_list = [0]
    num_atoms_last = 0

    for i in range(1, len(ag_list)):
        add_ag = ag_list[i].ag
        num_atoms_new = add_ag.numAtoms()

        # total things
        ag_out += add_ag

        # print("this offset:", i, num_atoms_last, len(bg_out))
        bg_list_offset = [x + offset for x in bg_list[i]]
        # different offset for connecting to the 5' end base
        if i > 1:
            offset1, offset2 = O3_ID_OFFSET_INT, O3_IDX_OFFSET_INT
        else:
            offset1, offset2 = O3_ID_OFFSET, O3_IDX_OFFSET

        # Create Phosphate Bond 3' to 5'
        base_from = offset__idx_list[-1] + offset1
        bg_list_offset[0][1] = base_from
        # Create Phosphate Bond 5' to 3'
        base_to = bg_list_offset[0][0]
        bg_out[-num_atoms_last + offset2][2] = base_to
        print("Connecting (%d, %d)" % (base_from, base_to))
        # print(-num_atoms_last + offset2, num_atoms_last, offset2, len(bg_out))

        # print("off last", num_atoms_last)
        bg_out += bg_list_offset


        offset__idx_list.append(offset)
        offset += num_atoms_new
        num_atoms_last = num_atoms_new
    # end for

    """
    install phosphate bond by 
    1. dropping the last Atom of the 5' side and
    2. the -OH atoms of the of the 3' side. 
    3. drop the last bond of the 5' side (the HCAP to the O) (base_bond_groups[i][-1])
    4. for the 3' side, replace the base_bond_groups[i+1][1] of the 3' 
        base_bond_groups[i+1][2][1] = offset = 
    """

    ag_out.setTitle(name)
    out_file = 'test_file.pdb'
    writePDB(out_file, ag_out)
    # print("bgOut len:", len(bg_out))
    writePDBConnect(out_file, bg_out)

    if create_psf:
        num_bonds = sum(len(x) for x in bg_out) - len(bg_out)
        bonds = []
        maxi = 0
        for item in bg_out:
            i = item[0]
            for j in item[1:]:
                if i < j: # so we don't double count
                    bonds.append((i-1,j-1))
                if i > maxi:
                    maxi = i
                if j > maxi:
                    maxi = j
        print(ag_out.numAtoms(), maxi)
        ag_out.setBonds(bonds)
        writePSF('test_file.psf', ag_out)
    return ag_out
# end def

if __name__ == "__main__":
    createStrand("ACGTACGTACG", None, create_psf=True)
    