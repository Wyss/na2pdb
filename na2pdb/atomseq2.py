import os.path as op
from collections import namedtuple
try:
    import na2pdb.matrix as matrix
    from na2pdb.cml import writeCML
    from na2pdb.mmcif import writeMMCIF
    from na2pdb.pdbfile import parsePDB, writePDB
    from na2pdb.parsebonds import parsePDBConnect, writePDBConnect
except:
    import matrix
    from cml import writeCML
    from mmcif import writeMMCIF
    from pdbfile import parsePDB, writePDB
    from parsebonds import parsePDBConnect, writePDBConnect
import math
import numpy as np


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
                    'T_single2.pdb'
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
DELTA_X_REV_OFFSET = 0.78 # unclear, but required

PRETWIST_Y_OFFSET = 2.29005198e-03 # for reverse strands NOT SURE how to use this
M_Y_OFF = matrix.makeTranslation(0, PRETWIST_Y_OFFSET, 0)

RADIUS = 10.175   # Approximate radius of each base in Angstroms, unused

THETA_REV_OFFSET = 0.18159025402704151

BASE_LUT = np.zeros(128, dtype=int)
BASE_LUT[b'A'[0]] = 0
BASE_LUT[b'C'[0]] = 1
BASE_LUT[b'G'[0]] = 2
BASE_LUT[b'T'[0]] = 3
BASE_LUT[b'a'[0]] = 0
BASE_LUT[b'c'[0]] = 1
BASE_LUT[b'g'[0]] = 2
BASE_LUT[b't'[0]] = 3

ROT = np.array([[ -1.0, -2.22044605e-16, -2.22044605e-16, 4.16],
                [ -6.17095915e-03, -9.83378029e-01, -1.82110602e-01, 1.89824119e-03],
                [  1.41609991e-04, -1.80761798e-01, 9.83150573e-01, 1.58762083e-03],
                [  0.0,   0.0,  5.55111512e-17,   1.0]])

def loadBasePDBs():
    base_atom_groups = []
    base_bonds = []
    for fn in dna_pdb_files:
        the_file = op.join(DATA_DIR, fn)
        ag = parsePDB(the_file)
        # print(ag.getDataLabels())
        base_atom_groups.append(ag)

        sub_bonds = parsePDBConnect(the_file)
        base_bonds.append(sub_bonds)
    return base_atom_groups, base_bonds

BASE_ATOM_GROUPS, BASE_BONDS = loadBasePDBs()

def create5PrimeAGs(base_atom_groups, base_bonds):
    """ Drop the trailing HCAP Hydrogen
    """
    five_p_atom_groups = []
    bonds = []
    for i, ag in enumerate(base_atom_groups):
        ag_copy = ag.copy()

        # ag_copy._n_atoms = ag_copy.numAtoms() - 1

        end_idx = ID_OFFSETS[i][1]

        ag_copy.df = ag_copy.df[:end_idx].reset_index(drop=True)
        new_coords = ag_copy._getCoords()[:end_idx]
        ag_copy._coords = None
        ag_copy.setCoords(new_coords)
        # print("coord check %d (%d, %d) %d" %(ag_copy._n_atoms, len(ag._coords[0]), len(ag_copy._coords[0]), end_idx))
        five_p_atom_groups.append(ag_copy)
        bonds.append([x.copy() for x in base_bonds[i][:end_idx]])
    return five_p_atom_groups, bonds
# end def
FIVE_P_ATOM_GROUPS, FIVE_P_BONDS = create5PrimeAGs(BASE_ATOM_GROUPS, BASE_BONDS)

def create3PrimeAGs(base_atom_groups, base_bonds):
    """ Drop the lead OH atoms
    """
    three_p_atom_groups = []
    bonds = []
    for i, ag in enumerate(base_atom_groups):
        ag_copy = ag.copy()

        # ag_copy._n_atoms = ag_copy.numAtoms() - 2

        start_idx = IDX_OFFSETS[i][0]

        ag_copy.df = ag_copy.df[start_idx:].reset_index(drop=True)
        
        new_coords = ag_copy._getCoords()[start_idx:]
        ag_copy._coords = None
        ag_copy.setCoords(new_coords)
        
        three_p_atom_groups.append(ag_copy)
        bonds.append([x - 2 for x in base_bonds[i][start_idx:]])
    return three_p_atom_groups, bonds
# end def
THREE_P_ATOM_GROUPS, THREE_P_BONDS = create3PrimeAGs(BASE_ATOM_GROUPS, BASE_BONDS)

def createINTAGs(base_atom_groups, base_bonds):
    """ Drop the lead OH atoms and trailing HCAP Hydrogen
    For internal bases
    """
    int_atom_groups = []
    bonds = []
    for i, ag in enumerate(base_atom_groups):
        ag_copy = ag.copy()
        # ag_copy._n_atoms = ag_copy.numAtoms() - 3

        start_idx = IDX_OFFSETS[i][0]
        end_idx = ID_OFFSETS[i][1]

        ag_copy.df = ag_copy.df[start_idx:end_idx].reset_index(drop=True)

        new_coords = ag_copy._getCoords()[start_idx:end_idx]
        ag_copy._coords = None
        ag_copy.setCoords(new_coords)

        int_atom_groups.append(ag_copy)
        bonds.append([x - 2 for x in base_bonds[i][start_idx:end_idx]])
    return int_atom_groups, bonds
# end def
INT_ATOM_GROUPS, INT_BONDS = createINTAGs(BASE_ATOM_GROUPS, BASE_BONDS)

class AtomicSequence(object):
    """
    1. construct all bases
    2. move bases to correct positions
    3. Apply twist.
    """
    def __init__(self, seq, name='', bases_per_turn=10.5, theta_offset=0.0):
        self.seq = seq
        self.bases_per_turn = bases_per_turn
        self.theta_offset = theta_offset
        self.atom_group = None

        # list of bonds
        self.bonds = None

        # list of virtual helix positions per bases where the  
        self.base_idxs = list(range(len(seq)))
        twist_per_segment = 2.*math.pi/self.bases_per_turn
        self.twists = [q*twist_per_segment + theta_offset for q in range(len(seq))]

        # list of index offsets for the atoms in the combined group
        self.start_idxs = [0]
        start_idxs = self.start_idxs

        self.reverse_queue = []
        self.transform_queue = []

        if isinstance(seq, str):
            seq = seq.encode('utf-8')
        ag_list = []
        bond_list = []
        lim = len(seq) - 1
        offset = 0
        for i, base in enumerate(seq):
            idx = BASE_LUT[base]
            if i == 0:
                ag = FIVE_P_ATOM_GROUPS[idx].copy()
                # copy a list a lists
                # FIVE_P_BONDS : [ [array(), array(),] ]
                bg = [x.copy() for x in FIVE_P_BONDS[idx]]
                offset += ag.numAtoms()
                start_idxs.append(offset)
            elif i < lim:
                ag = INT_ATOM_GROUPS[idx].copy()
                bg = [x.copy() for x in INT_BONDS[idx]]
                offset += ag.numAtoms()
                start_idxs.append(offset)
            else:
                ag = THREE_P_ATOM_GROUPS[idx].copy()
                bg = [x.copy() for x in THREE_P_BONDS[idx]]

            ag_list.append(AGBase(idx, ag))
            bond_list.append(bg)

        # for x in ag_list:
        #     print(x)

        """ Concatenate AtomGroups and Bonds apply offsets to bond indices
        reset offset variable will save a lookup into start_idxs
        """
        # prody.AtomGroup doesn't implement __radd__ so we can't do sum
        ag_out = ag_list[0].ag
        bonds_out = bond_list[0]
        offset = ag_out.numAtoms()
        num_atoms_last = 0

        for i in range(1, len(ag_list)):
            add_ag = ag_list[i].ag
            num_atoms_new = add_ag.numAtoms()
            # increment the residue sequence number to prevent atom name clashes
            add_ag.df['resSeq'] = 1 + i
            # total things
            ag_out.concat(add_ag)

            bond_list_offset = [x + offset for x in bond_list[i]]
            # different offset for connecting to the 5' end base
            if i > 1:
                offset1, offset2 = O3_ID_OFFSET_INT, O3_IDX_OFFSET_INT
            else:
                offset1, offset2 = O3_ID_OFFSET, O3_IDX_OFFSET

            # Create Phosphate Bond 3' to 5'
            base_from = start_idxs[i - 1] + offset1
            bond_list_offset[0][1] = base_from
            # Create Phosphate Bond 5' to 3'
            base_to = bond_list_offset[0][0]
            bonds_out[-num_atoms_last + offset2][2] = base_to
            # print("Connecting (%d, %d)" % (base_from, base_to))

            bonds_out += bond_list_offset
            offset += num_atoms_new
            num_atoms_last = num_atoms_new
        # end for

        """
        install phosphate bond by 
        1. dropping the last Atom of the 5' side and
        2. the -OH atoms of the of the 3' side. 
        3. drop the last bond of the 5' side (the HCAP to the O) (base_bonds[i][-1])
        4. for the 3' side, replace the base_bonds[i+1][1] of the 3' 
            base_bonds[i+1][2][1] = offset = 
        """
        ag_out.name = name
        self.atom_group = ag_out
        self.bonds = bonds_out
    # end def

    def concat(self, other):
        if not isinstance(other, AtomicSequence ):
            raise TypeError('unsupported operand type(s) for +: {0} and '
                            '{1}'.format(repr(type(self).__name__),
                                         repr(type(other).__name__)))

        atom_offset = self.atom_group.numAtoms()
        name = self.atom_group.getTitle()
        # 1. offset bonds
        obonds = other.bonds
        for i in range(len(other.bonds)):
            obonds[i] += atom_offset    # add scalar to array
        # add atom groups
        self.atom_group.concat(other.atom_group)
        self.atom_group.setTitle(name)
        self.bonds += other.bonds
        # don't bother with self.start_idxs, self.twists, etc for now
    # end def

    def transformBases(self, start, end, x, y, z, is_5to3):
        """
        assumes start < end
        start (int): the start index to transform
        end (int): the end_idx to transform (use -1 for the end) (non-inclusive)
        x (int): unit in bases from 0
        y, z (float): units in multiple of RADIUS from 0,0
        is_5to3 (bool): if the bases will be 5' to 3' in the x direction
            or 3' to 5'
        """
        old_coords = self.atom_group._getCoords()   # points to source array
        start_idxs = self.start_idxs
        twist_per_segment = 2.*math.pi/self.bases_per_turn
        theta0 = self.theta_offset

        start_idx = self.start_idxs[start]
        if end == -1 or end > len(self.seq) - 1:
            # print("doop", end)
            end_idx = len(old_coords)
        else:
            # print("goop", end)
            end_idx = self.start_idxs[end]
            
        if not is_5to3:
            # 1. Flip 180 degrees about Z to change direction
            m_rev = matrix.makeRotationZ(math.pi)
            self.reverse_queue.append((m_rev, start_idx, end_idx))
            # self.reverse_queue.append((M_Y_OFF, start_idx, end_idx))  

            # 2. Translate as required
            # print("translating", x + end - start )
            m = matrix.makeTranslation((x + end - start)*DELTA_X + DELTA_X_REV_OFFSET, 
                                        y*RADIUS, 
                                        z*RADIUS)
            if x + end - start - 1 < 1:
                print("tb", len(self.twists), x, end, start, x + end - start - 1)
                raise ValueError("poop")
            self.twists[start:end] = [(q*twist_per_segment + theta0 + THETA_REV_OFFSET)\
                                         for q in range(x + end - start - 1, x - 1, -1)]
        else:
            m = matrix.makeTranslation(x*DELTA_X, y*RADIUS, z*RADIUS)

            self.twists[start:end] = \
                            [(q*twist_per_segment + theta0) \
                                for q in range(x, x + end - start)]

        self.base_idxs[start:end] = list(range(0, end - start))

        # print(self.base_idxs)

        self.transform_queue.append((m, start_idx, end_idx))
    # end def

    def applyReverseQueue(self):
        new_coords = self.atom_group._getCoords()
        for item in self.reverse_queue:
            m, start, end = item
            new_coords[start:end] = matrix.applyTransform(new_coords[start:end], m)
        self.reverse_queue = []
    # end def

    def applyTransformQueue(self):
        new_coords = self.atom_group._getCoords()
        for item in self.transform_queue:
            m, start, end = item
            new_coords[start:end] = matrix.applyTransform(new_coords[start:end], m)
        self.transform_queue = []
    # end def

    def linearize(self):
        """  Using self.base_idxs, separate out bases relative to one another
        """
        new_coords = self.atom_group._getCoords()
        bidxs = self.base_idxs
        sidxs = self.start_idxs
        start = 0
        lim = len(sidxs) - 1
        for i in range(lim+1):
            if i < lim:
                next = sidxs[i+1]
            else:
                next = len(new_coords)
            base_idx = bidxs[i]
            m = matrix.makeTranslation(base_idx*DELTA_X, 0, 0)
            new_coords[start:next] = matrix.applyTransform(new_coords[start:next], m)
            start = next 
        # end for
        self.atom_group.setCoords(new_coords)
    # end def

    def applyTwist(self):
        """ Using self.base_idxs, twist bases relative to one another
        """
        new_coords = self.atom_group._getCoords()
        tidxs = self.twists
        sidxs = self.start_idxs
        # print("checker", len(tidxs), len(sidxs))
        assert(len(tidxs) == len(sidxs))
        start = 0
        lim = len(sidxs) - 1
        for i in range(lim+1):
            if i < lim:
                next = sidxs[i+1]
            else:
                next = len(new_coords)
            # print("twists", i, len(tidxs))
            theta = tidxs[i]
            m = matrix.makeRotationX(theta)
            new_coords[start:next] = matrix.applyTransform(new_coords[start:next], m)
            start = next
        # end for
        self.atom_group.setCoords(new_coords)
    # end def

    def toPDB(self, filename):
        ag_out = self.atom_group
        writePDB(filename, ag_out)
        bonds_out = self.bonds
        writePDBConnect(filename, bonds_out)
        writeCML(filename + '.cml', self)
        writeMMCIF(filename + '.cif', self)
    # end def

    def toPSF(self, filename):
        ag_out = self.atom_group
        bonds_out = self.bonds
        num_bonds = sum(len(x) for x in bonds_out) - len(bonds_out)
        bonds = []
        maxi = 0
        for item in bonds_out:
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
        writePSF(filename, ag_out)
# end class


def createStrand(seq, 
                name="strand", 
                bases_per_turn=10.5,
                theta_offset=0.0,
                create_psf=False):
    """
    """
    atom_sequence = AtomicSequence(seq, name=name, 
                                    bases_per_turn=bases_per_turn,
                                    theta_offset=theta_offset)

    len_strand = len(seq)
    half_len = len_strand // 2
    atom_sequence.transformBases(0, half_len, half_len, 0, 0, True)
    atom_sequence.transformBases(half_len, len_strand, half_len, 0, 0, False)
    # 1. Get base separation
    atom_sequence.linearize()
    # 2. do all rotations
    atom_sequence.applyReverseQueue()
    atom_sequence.applyTwist()
    # 3. move to position
    atom_sequence.applyTransformQueue()

    out_file = 'test_file.pdb'
    atom_sequence.toPDB(out_file)

    if create_psf:
        atom_sequence.toPSF('test_file.psf')
    return atom_sequence.atom_group
# end def

if __name__ == "__main__":
    # createStrand("ACGTACGTACG", None, create_psf=True)
    # a = R*b
    a = np.eye(4)
    a[:,0] = 4.16, -8.76, -1.609, 1    # P phosphate
    a[:,1] = 3.94, -10.139, -1.116, 1  # O1P
    a[:,2] = 3.37, -8.365, -2.796, 1   # O2P
    a[:,3] = 3.91, -7.718, -0.43, 1    # O5'

    b = np.eye(4)
    b[:,0] = 0., 8.91, 0, 1            # P phosphate
    b[:,1] = 0.22, 10.175, 0.734, 1    # O1P
    b[:,2] = 0.79, 8.733, -1.24, 1     # O2P
    b[:,3] = 0.25, 7.669, 0.971, 1     # O5'

    R = np.dot(a, np.linalg.inv(b))
    # print(R)
    # print(np.dot(R, b[:,0]))

    # print(np.dot(R, np.array([[-0.690, 7.424, 2.047, 1.00]]).T))
    # [[4.850,  -7.669,  0.674,  1.00]]
    createStrand("ACGTACGTACGTACGT", None)
    # createStrand("ATAT", None)



    
