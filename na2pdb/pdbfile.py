# -*- coding: utf-8 -*-
"""This module defines functions for parsing and writing `PDB files`_.

.. _PDB files: http://www.wwpdb.org/documentation/format32/v3.2.html

Derived from ProDy pdbfile.py but so as to use this version of AtomGroup

"""

from collections import defaultdict
import os.path

import numpy as np
STR_DTYPE = np.array(['a']).dtype.char  # 'S' for PY2K and 'U' for PY3K
ATOMIC_FIELDS = {    
    'name':     STR_DTYPE + '6',
    'altloc':   STR_DTYPE + '1',
    'anisou':   float,
    'chainid':  STR_DTYPE + '1',
    'element':  STR_DTYPE + '2',
    'occupancy':float,
    'resname':  STR_DTYPE + '6',
    'resnum':   int,
    'secondary':STR_DTYPE + '1',
    'segment':  STR_DTYPE + '6',
    'siguij':   float,
    'serial':   int,
    'beta':     float,
    'icode':    STR_DTYPE + '1',
    'type':     STR_DTYPE + '6',
    'charge':   STR_DTYPE + '2',
    'mass':     float,
    'radius':   float,
    'resindex': int,
    'chindex':  int, 
    'segindex': int,
    'fragindex':int,
    'numbonds': int
}

import numpy as np
from na2pdb.atomgroup import AtomGroup
import logging as LOGGER

__all__ = ['parsePDBStream', 'parsePDB',
           'writePDBStream', 'writePDB',]

class PDBParseError(Exception):
    pass


_parsePDBdoc = _parsePQRdoc + """
    :arg model: model index or None (read all models), e.g. ``model=10``
    :type model: int, list

    :arg altloc: if a location indicator is passed, such as ``'A'`` or ``'B'``,
         only indicated alternate locations will be parsed as the single
         coordinate set of the AtomGroup,  if *altloc* is set **True** all
         alternate locations will be parsed and each will be appended as a
         distinct coordinate set, default is ``"A"``
    :type altloc: str

    Note that this function does not evaluate ``CONECT`` records.

    """

_PDBSubsets = {'ca': 'ca', 'calpha': 'ca', 'bb': 'bb', 'backbone': 'bb'}

def parsePDB(pdb, **kwargs):
    """Return an :class:`.AtomGroup` and/or dictionary containing header data
    parsed from a PDB file.

    This function extends :func:`.parsePDBStream`.

    See :ref:`parsepdb` for a detailed usage example.

    :arg pdb: a PDB identifier or a filename
    """
    title, ext = os.path.splitext(os.path.split(pdb)[1])
    if ext == '.gz':
        title, ext = os.path.splitext(title)
    if len(title) == 7 and title.startswith('pdb'):
        title = title[3:]
    kwargs['title'] = title
    with open(pdb, 'r') as fd:
        result = parsePDBStream(fd, **kwargs)
    return result

parsePDB.__doc__ += _parsePDBdoc

def parsePDBStream(stream, **kwargs):
    """Return an :class:`.AtomGroup` and/or dictionary containing header data
    parsed from a stream of PDB lines.

    :arg stream: Anything that implements the method ``readlines``
        (e.g. :class:`file`, buffer, stdin)"""

    model = kwargs.get('model')
    header = kwargs.get('header', False)
    assert isinstance(header, bool), 'header must be a boolean'

    altloc = kwargs.get('altloc', 'A')
    if model is not None:
        if isinstance(model, int):
            if model < 0:
                raise ValueError('model must be greater than 0')
        else:
            raise TypeError('model must be an integer, {0} is invalid'
                            .format(str(model)))
    title_suffix = ''

    ag = None
    if 'ag' in kwargs:
        ag = kwargs['ag']
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.numCoordsets()
    elif model != 0:
        ag = AtomGroup(str(kwargs.get('title', 'Unknown')) + title_suffix)
        n_csets = 0

    if model != 0:
        try:
            lines = stream.readlines()
        except AttributeError as err:
            try:
                lines = stream.read().split('\n')
            except AttributeError:
                raise err
        if not len(lines):
            raise ValueError('empty PDB file or stream')

        _parsePDBLines(ag, lines, split, model, altloc)
        if ag.numAtoms() > 0:
            LOGGER.debug('{0} atoms and {1} coordinate set(s) were '
                          'parsed in %.2fs.'.format(ag.numAtoms(),
                           ag.numCoordsets() - n_csets))
        else:
            ag = None
            LOGGER.warning('Atomic data could not be parsed, please '
                        'check the input file.')
    if model != 0:
        return ag

parsePDBStream.__doc__ += _parsePDBdoc

def _parsePDBLines(atomgroup, lines, split, model, altloc_torf):
    """Return an AtomGroup. See also :func:`.parsePDBStream()`.

    :arg lines: PDB/PQR lines
    :arg split: starting index for coordinate data lines"""

    only_subset = False
=
    onlycoords = False
    n_atoms = atomgroup.numAtoms()
    if n_atoms > 0:
        asize = n_atoms
    else:
        # most PDB files contain less than 99999 atoms
        asize = min(len(lines) - split, 99999)
    addcoords = False
    if atomgroup.numCoordsets() > 0:
        addcoords = True
    alength = asize
    coordinates = np.zeros((asize, 3), dtype=float)
    atomnames = np.zeros(asize, dtype=ATOMIC_FIELDS['name'])
    resnames = np.zeros(asize, dtype=ATOMIC_FIELDS['resname'])
    resnums = np.zeros(asize, dtype=ATOMIC_FIELDS['resnum'])
    chainids = np.zeros(asize, dtype=ATOMIC_FIELDS['chainid'])
    hetero = np.zeros(asize, dtype=bool)
    termini = np.zeros(asize, dtype=bool)
    altlocs = np.zeros(asize, dtype=ATOMIC_FIELDS['altloc'])
    icodes = np.zeros(asize, dtype=ATOMIC_FIELDS['icode'])
    serials = np.zeros(asize, dtype=ATOMIC_FIELDS['serial'])

    segnames = np.zeros(asize, dtype=ATOMIC_FIELDS['segment'])
    elements = np.zeros(asize, dtype=ATOMIC_FIELDS['element'])
    charges = np.zeros(asize, dtype=ATOMIC_FIELDS['charge'])
    bfactors = np.zeros(asize, dtype=ATOMIC_FIELDS['beta'])
    occupancies = np.zeros(asize, dtype=ATOMIC_FIELDS['occupancy'])
    anisou = None
    siguij = None

    asize = 2000 # increase array length by this much when needed

    start = split
    stop = len(lines)
    nmodel = 0
    # if a specific model is requested, skip lines until that one
    if model is not None and model != 1:
        for i in range(split, len(lines)):
            if lines[i][:5] == 'MODEL':
                nmodel += 1
                if model == nmodel:
                    start = i+1
                    stop = len(lines)
                    break
        if nmodel != model:
            raise PDBParseError('model {0} is not found'.format(model))
    if isinstance(altloc_torf, str):
        if altloc_torf.strip() != 'A':
            LOGGER.info('Parsing alternate locations {0}.'
                        .format(altloc_torf))
            which_altlocs = ' ' + ''.join(altloc_torf.split())
        else:
            which_altlocs = ' A'
        altloc_torf = False
    else:
        which_altlocs = ' A'
        altloc_torf = True

    acount = 0
    altloc = defaultdict(list)
    i = start
    END = False
    while i < stop:
        line = lines[i]
        startswith = line[0:6]

        if startswith == 'ATOM  ' or startswith == 'HETATM':
            if only_subset:
                atomname = line[12:16].strip()
                resname = line[17:21].strip()
                if not (atomname in subset and resname in protein_resnames):
                    i += 1
                    continue
            else:
                atomname = line[12:16]
                resname = line[17:21]

            chid = line[21]

            alt = line[16]
            if alt not in which_altlocs:
                altloc[alt].append((line, i))
                i += 1
                continue
            try:
                coordinates[acount, 0] = line[30:38]
                coordinates[acount, 1] = line[38:46]
                coordinates[acount, 2] = line[46:54]
            except:
                if acount >= n_atoms > 0:
                    if nmodel ==0:
                        raise ValueError(format + 'file and AtomGroup ag must '
                                         'have same number of atoms')
                    LOGGER.warning('Discarding model {0}, which contains more '
                            'atoms than first model does.'.format(nmodel+1))
                    acount = 0
                    nmodel += 1
                    coordinates = np.zeros((n_atoms, 3), dtype=float)
                    while lines[i][:6] != 'ENDMDL':
                        i += 1
                else:
                    raise PDBParseError('invalid or missing coordinate(s) at '
                                         'line {0}.'.format(i+1))
            if onlycoords:
                acount += 1
                i += 1
                continue

            try:
                serials[acount] = line[6:11]
            except ValueError:
                try:
                    serials[acount] = int(line[6:11], 16)
                except ValueError:
                    LOGGER.warning('Failed to parse serial number in line {0}.'
                                .format(i))
                    serials[acount] = serials[acount-1]+1
            altlocs[acount] = alt
            atomnames[acount] = atomname
            resnames[acount] = resname
            chainids[acount] = chid
            resnums[acount] = line[22:26]#.split()[0])
            icodes[acount] = line[26]

            try:
                occupancies[acount] = line[54:60]
            except:
                LOGGER.warning('failed to parse occupancy at line {0}'
                            .format(i))
            try:
                bfactors[acount] = line[60:66]
            except:
                LOGGER.warning('failed to parse beta-factor at line {0}'
                            .format(i))
            hetero[acount] = startswith[0] == 'H'
            segnames[acount] = line[72:76]
            elements[acount] = line[76:78]
            charges[acount] = line[78:80]

            acount += 1
            if n_atoms == 0 and acount >= alength:
                # if arrays are short extend them with zeros
                alength += asize
                coordinates = np.concatenate(
                    (coordinates, np.zeros((asize, 3), float)))
                atomnames = np.concatenate((atomnames,
                    np.zeros(asize, ATOMIC_FIELDS['name'])))
                resnames = np.concatenate((resnames,
                    np.zeros(asize, ATOMIC_FIELDS['resname'])))
                resnums = np.concatenate((resnums,
                    np.zeros(asize, ATOMIC_FIELDS['resnum'])))
                chainids = np.concatenate((chainids,
                    np.zeros(asize, ATOMIC_FIELDS['chainid'])))
                hetero = np.concatenate((hetero, np.zeros(asize, bool)))
                termini = np.concatenate((termini, np.zeros(asize, bool)))
                altlocs = np.concatenate((altlocs,
                    np.zeros(asize, ATOMIC_FIELDS['altloc'])))
                icodes = np.concatenate((icodes,
                    np.zeros(asize, ATOMIC_FIELDS['icode'])))
                serials = np.concatenate((serials,
                    np.zeros(asize, ATOMIC_FIELDS['serial'])))
                bfactors = np.concatenate((bfactors,
                    np.zeros(asize, ATOMIC_FIELDS['beta'])))
                occupancies = np.concatenate((occupancies,
                    np.zeros(asize, ATOMIC_FIELDS['occupancy'])))
                segnames = np.concatenate((segnames,
                    np.zeros(asize, ATOMIC_FIELDS['segment'])))
                elements = np.concatenate((elements,
                    np.zeros(asize, ATOMIC_FIELDS['element'])))
                charges = np.concatenate((charges,
                    np.zeros(asize, ATOMIC_FIELDS['charge'])))
                if anisou is not None:
                    anisou = np.concatenate((anisou, np.zeros((asize, 6),
                        ATOMIC_FIELDS['anisou'])))
                if siguij is not None:
                    siguij = np.concatenate((siguij, np.zeros((asize, 6),
                        ATOMIC_FIELDS['siguij'])))
        #elif startswith == 'END   ' or startswith == 'CONECT':
        #    i += 1
        #    break
        elif not onlycoords and (startswith == 'TER   ' or
            startswith.strip() == 'TER'):
            termini[acount - 1] = True
        elif startswith == 'ENDMDL' or startswith[:3] == 'END':
            if acount == 0:
                # If there is no atom record between ENDMDL & END skip to next
                i += 1
                continue
            if model is not None:
                i += 1
                break
            diff = stop - i - 1
            if diff < acount:
                END = True
            if onlycoords:
                if acount < n_atoms:
                    LOGGER.warning('Discarding model {0}, which contains '
                                '{1} fewer atoms than the first model '
                                'does.'.format(nmodel+1, n_atoms-acount))
                else:
                    coordsets[nmodel] = coordinates
                    nmodel += 1
                acount = 0
                if not END:
                    coordinates = coordsets[nmodel]
            else:
                if acount != n_atoms > 0:
                    raise ValueError('PDB file and AtomGroup ag must have '
                                    'same number of atoms')
                # this is where to decide if more coordsets should be expected
                if END:
                    coordinates.resize((acount, 3))
                    if addcoords:
                        atomgroup.addCoordset(coordinates)
                    else:
                        atomgroup._setCoords(coordinates)
                else:
                    coordsets = np.zeros((diff/acount+1, acount, 3))
                    coordsets[0] = coordinates[:acount]
                    onlycoords = True
                atomnames.resize(acount)
                resnames.resize(acount)
                resnums.resize(acount)
                chainids.resize(acount)
                hetero.resize(acount)
                termini.resize(acount)
                altlocs.resize(acount)
                icodes.resize(acount)
                serials.resize(acount)
                if not only_subset:
                    atomnames = np.char.strip(atomnames)
                    resnames = np.char.strip(resnames)

                bfactors.resize(acount)
                occupancies.resize(acount)
                segnames.resize(acount)
                elements.resize(acount)
                segnames = np.char.strip(segnames)

                if anisou is not None:
                    anisou.resize((acount, 6))
                    atomgroup.setAnisous(anisou / 10000)
                if siguij is not None:
                    siguij.resize((acount, 6))
                    atomgroup.setAnistds(siguij / 10000)

                atomgroup.setDataFrame(atomnames, altlocs, resnames, 
                                        chainids, resnums, np.char.strip(icodes),
                                        occupancies, bfactors, np.char.strip(elements),
                                        np.char.strip(charges))

                nmodel += 1
                n_atoms = acount
                acount = 0
                coordinates = np.zeros((n_atoms, 3), dtype=float)
                if altloc and altloc_torf:
                    _evalAltlocs(atomgroup, altloc, chainids, resnums,
                                 resnames, atomnames)
                    altloc = defaultdict(list)
                
                if END:
                    break
        elif startswith == 'ANISOU':
            if anisou is None:
                anisou = True
                anisou = np.zeros((alength, 6),
                    dtype=ATOMIC_FIELDS['anisou'])
            try:
                index = acount - 1
                anisou[index, 0] = line[28:35]
                anisou[index, 1] = line[35:42]
                anisou[index, 2] = line[43:49]
                anisou[index, 3] = line[49:56]
                anisou[index, 4] = line[56:63]
                anisou[index, 5] = line[63:70]
            except:
                LOGGER.warning('failed to parse anisotropic temperature '
                    'factors at line {0}'.format(i))
        elif startswith =='SIGUIJ':
            if siguij is None:
                siguij = np.zeros((alength, 6),
                    dtype=ATOMIC_FIELDS['siguij'])
            try:
                index = acount - 1
                siguij[index, 0] = line[28:35]
                siguij[index, 1] = line[35:42]
                siguij[index, 2] = line[43:49]
                siguij[index, 3] = line[49:56]
                siguij[index, 4] = line[56:63]
                siguij[index, 5] = line[63:70]
            except:
                LOGGER.warning('failed to parse standard deviations of '
                    'anisotropic temperature factors at line {0}'.format(i))
        elif startswith =='SIGATM':
            pass
        i += 1
    if onlycoords:
        if acount == atomgroup.numAtoms():
            coordsets[nmodel] = coordinates
            nmodel += 1
        del coordinates
        coordsets.resize((nmodel, atomgroup.numAtoms(), 3))
        if addcoords:
            atomgroup.addCoordset(coordsets)
        else:
            atomgroup._setCoords(coordsets)
    elif not END:
        # this means last line was an ATOM line, so atomgroup is not decorated
        coordinates.resize((acount, 3))
        if addcoords:
            atomgroup.addCoordset(coordinates)
        else:
            atomgroup._setCoords(coordinates)
        atomnames.resize(acount)
        resnames.resize(acount)
        resnums.resize(acount)
        chainids.resize(acount)
        hetero.resize(acount)
        termini.resize(acount)
        altlocs.resize(acount)
        icodes.resize(acount)
        serials.resize(acount)
        if not only_subset:
            atomnames = np.char.strip(atomnames)
            resnames = np.char.strip(resnames)

        if anisou is not None:
            anisou.resize((acount, 6))
            # atomgroup.setAnisous(anisou / 10000)
        if siguij is not None:
            siguij.resize((acount, 6))
            # atomgroup.setAnistds(siguij / 10000)
        bfactors.resize(acount)
        occupancies.resize(acount)
        segnames.resize(acount)
        segnames = np.char.strip(segnames)
        elements.resize(acount)

        atomgroup.setDataFrame(atomnames, altlocs, resnames, 
                        chainids, resnums, np.char.strip(icodes),
                        occupancies, bfactors, np.char.strip(elements),
                        np.char.strip(charges))

    if altloc and altloc_torf:
        _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames)

    return atomgroup

def _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames):
    altloc_keys = list(altloc)
    altloc_keys.sort()
    indices = {}
    for key in altloc_keys:
        xyz = atomgroup.getCoords()
        success = 0
        lines = altloc[key]
        for line, i in lines:
            aan = line[12:16].strip()
            arn = line[17:21].strip()
            ach = line[21]
            ari = int(line[22:26].split()[0])
            rn, ids, ans = indices.get((ach, ari), (None, None, None))
            if ids is None:
                ids = indices.get(ach, None)
                if ids is None:
                    ids = (chainids == ach).nonzero()[0]
                    indices[ach] = ids
                ids = ids[resnums[ids] == ari]
                if len(ids) == 0:
                    LOGGER.warning("failed to parse altloc {0} at line {1}, "
                                "residue not present for altloc 'A'".format(
                                repr(key), i+1))
                    continue
                rn = resnames[ids[0]]
                ans = atomnames[ids]
                indices[(ach, ari)] = (rn, ids, ans)
            if rn != arn:
                LOGGER.warning("failed to parse altloc {0} at line {1}, "
                            "residue name mismatch (expected {2}, "
                            "parsed {3})".format(repr(key), i+1, repr(rn),
                                                   repr(arn)))
                continue
            index = ids[(ans == aan).nonzero()[0]]
            if len(index) != 1:
                LOGGER.warning("failed to parse altloc {0} at line {1}, atom"
                            " {2} not found in the residue"
                            .format(repr(key), i+1, repr(aan)))
                continue
            try:
                xyz[index[0], 0] = float(line[30:38])
                xyz[index[0], 1] = float(line[38:46])
                xyz[index[0], 2] = float(line[46:54])
            except:
                LOGGER.warning('failed to parse altloc {0} at line {1}, could'
                            ' not read coordinates'.format(repr(key), i+1))
                continue
            success += 1
        LOGGER.info('{0} out of {1} altloc {2} lines were parsed.'
                    .format(success, len(lines), repr(key)))
        if success > 0:
            LOGGER.info('Altloc {0} is appended as a coordinate set to the '
                        'atom group.'.format(repr(key), atomgroup.getTitle()))
            atomgroup.addCoordset(xyz)

# PDBLINE = ('{0:6s}{1:5d} {2:4s}{3:1s}'
#            '{4:4s}{5:1s}{6:4d}{7:1s}   '
#            '{8:8.3f}{9:8.3f}{10:8.3f}'
#            '{11:6.2f}{12:6.2f}      '
#            '{13:4s}{14:2s}\n')

PDBLINE_LT100K = ('%-6s%5d %-4s%1s%-4s%1s%4d%1s   '
                  '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                  '    %2s%2s\n')

PDBLINE_GE100K = ('%-6s%5x %-4s%1s%-4s%1s%4d%1s   '
                  '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                  '    %2s%2s\n')


_writePDBdoc = """

    :arg atomgroup: an object with atom and coordinate data

    """

def writePDBStream(stream, atomgroup, **kwargs):
    """Write *atoms* in PDB format to a *stream*.

    :arg stream: anything that implements a :meth:`write` method (e.g. file,
        buffer, stdout)"""

    remark = atomgroup.name
    df = atomgroup.df
    coordsets = atomgroup.coords

    n_atoms = atomgroup.numAtoms()
    occupancies = df['occupancy']
    bfactors = df['tempFactor']
    atomnames = df['names'].copy()
    
    for i, an in enumerate(atomnames):
        if len(an) < 4:
            atomnames[i] = ' ' + an

    altlocs = df['altLoc']

    resnames = df['resName']

    chainids = df['chainID']

    resnums = df['resSeq']

    icodes = df['iCode']

    hetero = ['ATOM'] * n_atoms
    # heteroflags = atoms._getFlags('hetatm')
    # if heteroflags is None:
    #     heteroflags = atoms._getFlags('hetero')
    # if heteroflags is not None:
    #     hetero = np.array(hetero, ATOMIC_FIELDS['hetero'])
    #     hetero[heteroflags] = 'HETATM'

    elements = df['elements']
    # elements = np.char.rjust(elements, 2)

    # segments = atoms._getSegnames()
    # if segments is None:
    #     segments = np.zeros(n_atoms, ATOMIC_FIELDS['segments'])

    charges = df['charges']

    stream.write('REMARK {0}\n'.format(remark))

    multi = len(coordsets) > 1
    write = stream.write
    for m, coords in enumerate(coordsets):
        pdbline = PDBLINE_LT100K
        if multi:
            write('MODEL{0:9d}\n'.format(m+1))
        for i, xyz in enumerate(coords):
            if i == 99999:
                pdbline = PDBLINE_GE100K
            write(pdbline % (hetero[i], i+1,
                         atomnames[i], altlocs[i],
                         resnames[i], chainids[i], resnums[i],
                         icodes[i],
                         xyz[0], xyz[1], xyz[2],
                         occupancies[i], bfactors[i],
                         elements[i], charges[i]))
        if multi:
            write('ENDMDL\n')
            altlocs = np.zeros(n_atoms, ATOMIC_FIELDS['altloc'])

writePDBStream.__doc__ += _writePDBdoc

def writePDB(filename, atomgroup):
    """Write *atoms* in PDB format to a file with name *filename* and return
    *filename*.  If *filename* ends with :file:`.gz`, a compressed file will
    be written."""
    with open(filename, 'w') as fd:
        writePDBStream(fd)
    return filename

writePDB.__doc__ += _writePDBdoc + """
    :arg autoext: when not present, append extension :file:`.pdb` to *filename*
"""
