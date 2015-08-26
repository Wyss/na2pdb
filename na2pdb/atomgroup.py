import pandas as pd
import numpy as np
import copy

PDB_COLS = [ 'name',
            'altLoc',
            'resName',
            'chainID',
            'resSeq',
            'iCode',
            'occupancy',
            'tempFactor',
            'element',
            'charge']

class AtomGroup(object):
    def __init__(self, name=''):
        """

        """
        self.name = name
        self.df = None
        self.coords = [[]]
        self.trajectories = []
        self.bonds = None
        self._flags = {}
    # end def

    def setDataFrame(self, *args):
        """ expects the numyp as
        """
        if len(args) != len(PDB_COLS):
            raise ValueError("needs more args")
        data_dict = dict(zip(PDB_COLS, args))
        # for key, val in data_dict.items():
        #     print(key, len(val))
        self.df = pd.DataFrame(data_dict)
    # end def

    def addCoordset(self, coords):
        self.coords.append(coords)
    # end def

    def getCoordset(self, coords, index=0):
        return self.coords[index]
    # end def

    def _getCoords(self):
        return self.coords[0]
    # end def

    def setCoords(self, coords):
        self.coords[0] = coords
    # end def

    def setCoordset(self, coords, index=0):
        self.coords[index] = coords
    # end def

    def getElements(self):
        return self.df['element']
    # end def

    def numCoordsets(self):
        return len(self.coords)
    # end def

    def copy(self):
        copy_ag = AtomGroup(self.name)
        copy_ag.df = self.df.copy()
        copy_ag.coords = copy.deepcopy(self.coords)
        copy_ag.trajectories = copy.deepcopy(self.trajectories)
        copy_ag.bonds = copy.deepcopy(self.bonds)
        copy_ag._flags = copy.deepcopy(self._flags)
        return copy_ag
    # end def

    def concat(self, other):
        """ combine two atom groups while int
        """
        if not isinstance(other, AtomGroup):
            raise TypeError('unsupported operand type(s) for +: {0} and '
                            '{1}'.format(repr(type(self).__name__),
                                         repr(type(other).__name__)))
        offset = self.numAtoms()

        self.df = self.df.append(other.df, ignore_index=True)

        if len(self.coords) != len(other.coords):
            raise ValueError("need to have same number of coordinate sets")
        for i in range(len(self.coords)):
            self.coords[i] = np.vstack((self.coords[i], other.coords[i]))

        # only concatinate trajectories if frame counts match
        frame_count = self.frameCount()
        if frame_count > 0 and frame_count == other.frameCount():
            for i in range(frame_count):
                self.trajectories[i] = np.vstack((self.trajectories[i], 
                                                other.trajectories[i]))
            # end for
        # end if

        if self.bonds is not None:
            self.bonds = np.concatenate((self.bonds, other.bonds + other), axis=0)

    # end def

    def numAtoms(self):
        return len(self.coords[0])
    # end def

    def frameCount(self):
        return len(self.trajectories)
    # end def

    def setBonds(self, bonds):
        """Set covalent bonds between atoms.  *bonds* must be a list or an
        array of pairs of indices.  All bonds must be set at once.  Bonding
        information can be used to make atom selections, e.g. ``"bonded to
        index 1"``.  See :mod:`.select` module documentation for details.
        Also, a data array with number of bonds will be generated and stored
        with label *numbonds*.  This can be used in atom selections, e.g.
        ``'numbonds 0'`` can be used to select ions in a system."""

        if isinstance(bonds, list):
            bonds = np.array(bonds, int)
        if bonds.ndim != 2:
            raise ValueError('bonds.ndim must be 2')
        if bonds.shape[1] != 2:
            raise ValueError('bonds.shape must be (n_bonds, 2)')
        if bonds.min() < 0:
            raise ValueError('negative atom indices are not valid')
        n_atoms = self.numAtoms()
        if bonds.max() >= n_atoms:
            raise ValueError('atom indices are out of range')
        bonds.sort(1)
        bonds = bonds[bonds[:, 1].argsort(), ]
        bonds = bonds[bonds[:, 0].argsort(), ]
        self.bonds = bonds
    # end def

    def numBonds(self):
        """Return number of bonds.  Use :meth:`setBonds` for setting bonds."""

        if self._bonds is not None:
            return self._bonds.shape[0]
        return 0
    # end def

    def getFlags(self, key):
        """Return a copy of atom flags for given *key*, or **None** when
        flags for *key* is not set."""

        flags = self._getFlags(key)
        if flags is not None:
            return flags.copy()

    def _getFlags(self, key):
        """Return atom flag values for given *key*, or **None** when
        flags for *key* is not set."""
        return self._flags[key]

    def setFlags(self, key, flags):
        """Set atom *flags* for *key*."""

        try:
            ndim, dtype = flags.ndim, flags.dtype
        except AttributeError:
            flags = np.array(flags)
            ndim, dtype = flags.ndim, flags.dtype
        if ndim != 1:
            raise ValueError('flags.ndim must be 1')
        if dtype != bool:
            raise ValueError('flags.dtype must be bool')
        if len(flags) != self.numAtoms():
            raise ValueError('len(flags) must be equal to number of atoms')
        self._setFlags(key, flags)

    def _setFlags(self, key, flags):
        """Set atom flags."""
        self._flags[key] = flags

    def delFlags(self, key):
        """Return flags associated with *key* and remove from the instance.
        If flags associated with *key* is not found, return **None**."""
        return self._flags.pop(key, None)
# end class