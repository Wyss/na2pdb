
from prody.proteins import parsePDB, writePDB
import os.path as op
from collections import namedtuple
import matrix
import numpy as np
"""
# strands axii assumed in the X direction

1. load bases pdbs
2. construct an list of prody.AtomGroups for the strand
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
IDX_OFFSETS = [ (3, 34),
                (3, 32),
                (3, 35),
                (3, 34)
            ]
# convert to zero-based indexing
IDX_OFFSETS = [(x - 1, y - 1) for x, y in IDX_OFFSETS]

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
    for fn in dna_pdb_files:
        ag = parsePDB(op.join(DATA_DIR, fn))
        base_atom_groups.append(ag)
    return base_atom_groups

base_atom_groups = loadBasePDBs()

def createStrand(seq, origin, name="strand", is_fwd=False, theta_offset=0.0):
    """
    """
    global base_atom_groups
    if isinstance(seq, str):
        seq = seq.encode('utf-8')
    ag_list = []
    for i, base in enumerate(seq):
        idx = BASE_LUT[base]
        print(chr(base), idx)
        ag = base_atom_groups[idx].copy()
        m = matrix.makeTranslation(i*DELTA_X, 0, 0)
        new_coords = matrix.applyTransform(ag._getCoords(), m) 
        ag.setCoords(new_coords)
        ag_list.append(AGBase(idx, ag))

    for x in ag_list:
        print(x)

    # prody.AtomGroup doesn't impplement __radd__ so we can't do sum
    ag_out = ag_list[0].ag
    for i in range(1, len(ag_list)):
        ag_out += ag_list[i].ag

    ag_out.setTitle(name)
    writePDB('test_file.pdb', ag_out)
    return ag_out
# end def


    
