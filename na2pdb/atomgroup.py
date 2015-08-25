import pandas as pd
import numpy as np

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
        self.coords = []
        self.trajectories = None
        self.bonds = None
    # end def

    def setDataFrame(self, *args):
        """ expects the numyp as
        """
        if len(args) != len(PDB_COLS):
            raise ValueError("needs more args")
        Series = 
        data_dict = dict(zip(PDB_COLS, args))
        self.df = pd.DataFrame(data_dict)
    # end def

    def addCoordset(self, coords):
        self.coords.appends(coords)
    # end def

    def getCoordset(self, coords, index=0):
        return self.coords[index]
    # end def

    def setCoordset(self, coords, index=0):
        self.coords[index] = coords
    # end def

    def concat(self, other):
        """ combine two atom groups while int
        """
        if not isinstance(other, AtomGroup):
            raise TypeError('unsupported operand type(s) for +: {0} and '
                            '{1}'.format(repr(type(self).__name__),
                                         repr(type(other).__name__)))
        offset = self.numAtoms()

        self.df.append(other.df, ignore_index=True)

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

        self.bonds = np.concatenate((self.bonds, other.bonds + other), axis=0)

    # end def

    """ getters and setters based on the PDB atom format
    """

    def numAtoms(self):
        return len(self.coords)
    # end def

    def frameCount(self):
        return len(self.trajectories)
    # end def

    def setNames(self):

    # end def

    def getNames(self):

    # end def

    def setAltlocs(self):

    # end def

    def getAltLocs(self):

    # end def

    def getResName(self):

    # end def

    def setResName(self):

    # end def

    def getChainID(self):

    # end def

    def setChainID(self):

    # end def

    def getResSeqs(self):

    # end def

    def setResSeqs(self):

    # end def

    def getCoords(self):

    # end def

    def setCoords(self):

    # end def

    def getOccupancies(self):

    # end def

    def setOccupancies(self):

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

    def __add__(self, other):



        new = AtomGroup(self._title + ' + ' + other._title)
        if self._n_csets:
            if self._n_csets == other._n_csets:
                new.setCoords(np.concatenate((self._coords, other._coords), 1))
                if self._n_csets > 1:
                    LOGGER.info('All {0} coordinate sets are copied to '
                                '{1}.'.format(self._n_csets, new.getTitle()))
            else:
                new.setCoords(np.concatenate((self._getCoords(),
                                              other._getCoords())))
                LOGGER.info('Active coordinate sets are copied to {0}.'
                            .format(new.getTitle()))
        elif other._n_csets:
            LOGGER.warn('No coordinate sets are copied to {0}'
                        .format(new.getTitle()))

        for key in set(list(self._data) + list(other._data)):
            if key in ATOMIC_FIELDS and ATOMIC_FIELDS[key].readonly:
                continue
            this = self._data.get(key)
            that = other._data.get(key)
            if this is not None or that is not None:
                if this is None:
                    this = np.zeros(that.shape, that.dtype)
                if that is None:
                    that = np.zeros(this.shape, this.dtype)
                new._data[key] = np.concatenate((this, that))

        if self._bonds is not None and other._bonds is not None:
            new.setBonds(np.concatenate([self._bonds,
                                         other._bonds + self._n_atoms]))
        elif self._bonds is not None:
            new.setBonds(self._bonds.copy())
        elif other._bonds is not None:
            new.setBonds(other._bonds + self._n_atoms)

        return new
# end class