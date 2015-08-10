import numpy as np

"""see http://www.bmsc.washington.edu/CrystaLinks/man/pdb/guide2.2_frame.html
We don't care about the details really other than the fact we'll need to add 
the phosphate bond and get rid  of a hydrogen bond and an OH bond per base
"""
def parsePDBConnect(fn):
    """ return a numpy array of the bonds so we can add an offset scalar to 
    its 
    """
    bonds = []
    with open(fn, 'r') as fd:
        for line in fd:
            startswith = line[0:6]
            if startswith == 'CONECT':
                bonds_per_atom = [int(x) for x in line[7:].split()]
                if len(bonds_per_atom) > 1:
                    bonds.append(np.array(bonds_per_atom, dtype=int))
    print("num_bond lines:", len(bonds))
    return bonds
# end def

def writePDBConnect(fn, bonds):
    """ append to a ProDy written pdb
    ProDy doesn't write an end
    """
    with open(fn, "a") as fd:
        for bonds_per_atom in bonds:
            num_bonds = len(bonds_per_atom)
            outstring = "CONECT"
            outstring += ''.join(["%5d" for x in range(num_bonds)])
            outstring += '\n'
            fd.write(outstring % tuple(bonds_per_atom))
        fd.write('END\n')
# end def