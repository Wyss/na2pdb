def writeCML(filename, aseq):
    bonds_out = aseq.bonds
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
    with open(filename, 'w') as fd:
        fd.write("<molecule>\n")
        fd.write(" <atomArray>\n")
        ag = aseq.atom_group
        coords = ag._getCoords()
        elements = ag.getElements()
        atom_str = "  <atom id=\"a%d\" elementType=\"%s\" x3=\"%f\" y3=\"%f\" z3=\"%f\"/>\n"
        for i in range(ag.numAtoms()):
            element = elements[i]
            x, y, z = coords[i]
            fd.write(atom_str % (i+1, element, x, y, z))
        fd.write(" </atomArray>\n")
        fd.write(" <bondArray>\n")
        bond_str = "<bond atomRefs2=\"a%d a%d\" order=\"1\"/>"
        for bond in bonds:
            fd.write(bond_str % (bond[0], bond[1]))
        fd.write(" </bondArray>\n")
        fd.write("</molecule>\n")
# end def