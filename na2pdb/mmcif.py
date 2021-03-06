"""
This is based upon requirements for successfully parsing
mmCIF files in Chimera and that the struct_conn_type_id 
for covalent bonds is `covale` as per:
 
http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_struct_conn_type.id.html

and not just a `c`
"""
ROOT_STR = \
"""
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.type_symbol
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.occupancy
_atom_site.pdbx_auth_alt_id
_atom_site.B_iso_or_equiv
"""

GROUP_STR = "ATOM"
ID_STR = " %d"
AUTH_ATOM_ID_STR = " %s"
TYPE_SYMBOL_STR = " %s"
AUTH_COMP_ID_STR = " %s"
AUTH_ASYM_ID_STR = " %s"
AUTH_SEQ_ID_STR = " %d"
POSITION_STR = " %f %f %f"
PDBX_PDB_MODEL_NUM_STR = " %d"
OCCUPANCY_STR = " %f"
PDBX_AUTH_ALT_ID_STR = " %s"
B_ISO_OR_EQUIV_STR = " %f"



ATOM_STR = GROUP_STR + ID_STR + AUTH_ATOM_ID_STR + \
            TYPE_SYMBOL_STR + AUTH_COMP_ID_STR + AUTH_ASYM_ID_STR +\
            AUTH_SEQ_ID_STR +\
            POSITION_STR + PDBX_PDB_MODEL_NUM_STR + OCCUPANCY_STR +\
            PDBX_AUTH_ALT_ID_STR + B_ISO_OR_EQUIV_STR + '\n'

STRUCT_CONN_HEADER =\
"""
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id 
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id 
_struct_conn.ptnr2_label_asym_id 
_struct_conn.ptnr2_label_comp_id 
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
"""
BOND_STR = "%d covale %s %s %d %s %s %s %d %s\n"

"""
1  c C DT 3  N3 D A 27 N1
2  c C DT 3  O4 D A 27 N6
3  c C DA 4  N1 D U 26 N3
4  c C DA 4  N6 D U 26 O4

http://www.cgl.ucsf.edu/chimera/data/mmcif-oct2013/mmcif.html
mandatory columns: chain id, residue name, residue number and atom name of joined atoms

"""

def writeMMCIF(filename, atomic_sequence):
    bonds_out = atomic_sequence.bonds
    num_connected_atoms = len(bonds_out)
    num_bonds = sum(len(x) for x in bonds_out) - num_connected_atoms
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
    with open(filename, 'w') as fd:
        ag = atomic_sequence.atom_group
        fd.write("data_%s" % ag.name)
        fd.write(ROOT_STR)
        coords = ag._getCoords()
        elements = ag.getElements()
        resnames = ag.df['resName']
        names = ag.df['name']
        chainids = ag.df['chainID']
        seqids = ag.df['resSeq']
        occupancies = ag.df['occupancy']
        temp_factors = ag.df['tempFactor']
        for i in range(ag.numAtoms()):
            element = elements[i]
            x, y, z = coords[i]
            name = names[i]
            resname = resnames[i]
            cid = chainids[i]
            seqid = seqids[i]
            model_num = 1
            fd.write(ATOM_STR % (i+1, name, element, resname, cid, seqid, 
                                x, y, z, model_num, occupancies[i], 
                                '.', temp_factors[i]))
        fd.write(STRUCT_CONN_HEADER)
        for i, bond in enumerate(bonds):
            b0, b1 = bond
            chainid0 = chainids[b0]
            rn0 = resnames[b0]
            seqid0 = seqids[b0]
            name0 = names[b0]

            chainid1 = chainids[b1]
            rn1 = resnames[b1]
            seqid1 = seqids[b1]
            name1 = names[b1]

            fd.write(BOND_STR % (i+1, chainid0, rn0, seqid0, name0, chainid1, rn1, seqid1, name1))
        # end for
# end def

# ATOM       1  N    N  GLN  A   39   24.690  -27.754   24.275  1  1.000  .  60.760
# ATOM       2  CA   C  GLN  A   39   23.581  -26.768   24.416  1  1.000  .  60.980
# ATOM       3  C    C  GLN  A   39   23.990  -25.379   23.905  1  1.000  .  59.980
# ATOM       4  O    O  GLN  A   39   25.070  -25.209   23.330  1  1.000  .  60.250
# ATOM       5  CB   C  GLN  A   39   23.136  -26.685   25.878  1  1.000  .  60.690
# ATOM       6  N    N  VAL  A   40   23.115  -24.395   24.122  1  1.000  .  59.580
# ATOM       7  CA   C  VAL  A   40   23.342  -23.010   23.690  1  1.000  .  57.260
# ATOM       8  C    C  VAL  A   40   24.000  -22.152   24.778  1  1.000  .  56.000
# ATOM       9  O    O  VAL  A   40   23.992  -20.920   24.692  1  1.000  .  55.530
# ATOM      10  CB   C  VAL  A   40   22.015  -22.337   23.275  1  1.000  .  57.320
# ATOM      11  N    N  ALA  A   41   24.560  -22.804   25.797  1  1.000  .  54.570