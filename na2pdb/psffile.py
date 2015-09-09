
"""
Derived from ProDy's trajectory.psffile
"""
PSFLINE = ('%8d %-4s %-4d %-4s %-4s %-4s %10.6f %13.4f %11d\n')

def writePSF(filename, atomseq):
    """Write atoms in X-PLOR format PSF file with name *filename* and return
    *filename*.  This function will write available atom and bond information
    only."""
    atomgroup = atomseq.atom_group
    bonds_out = atomseq.bonds
    num_bonds = sum(len(x) for x in bonds_out) - len(bonds_out)
    bonds = []
    maxi = 0
    for item in bonds_out:
        i = item[0]
        for j in item[1:]:
            if i < j: # so we don't double count
                bonds.append((i,j))
            if i > maxi:
                maxi = i
            if j > maxi:
                maxi = j

    try:
        n_atoms = atomgroup.numAtoms()
        df = atomgroup.df
        segments = chainids = df['chainID']
        elements = df['element']   # aka
        resnums = df['resSeq']
        resnames = df['resName']
        atomnames = df['name']
        atomtypes = None
        charges = None
        masses = None

    except AttributeError:
        raise TypeError('atoms must be an atomseq')

    if segments is None:
        segments = ['UNK'] * n_atoms

    if atomtypes is None:
        atomtypes = zeros(n_atoms, array(['a']).dtype.char + '1')

    if charges is None:
        charges = zeros(n_atoms, float)
    if masses is None:
        masses = zeros(n_atoms, int)

    long_fields = array([len(tp) for tp in types]).max() > 4

    with open(filename, 'w') as out:
        write = out.write
        write('PSF{0}\n'.format( ' NAMD' if long_fields else ''))
        write('\n')
        write('{0:8d} !NTITLE\n'.format(1))
        write(' REMARKS {0}\n'.format(atomgroup.name))
        write('\n')
        write('{0:8d} !NATOM\n'.format(n_atoms))

        for i in range(n_atoms):
            write(PSFLINE % (i + 1, segments[i], resnums[i],
                            resnames[i], atomnames[i], atomtypes[i], 
                            charges[i], masses[i], 0) )
        
        write('\n')
        write('{0:8d} !NBOND: bonds\n'.format(len(bonds)))
        for i, bond in enumerate(bonds):
            write('%8s%8s' % (bond[0], bond[1]))
            if i % 4 == 3:
                write('\n')
        if i % 4 != 3:
            write('\n')
    # end with
    return filename
# end def