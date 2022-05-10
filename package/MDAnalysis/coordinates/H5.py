# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
H5 trajectories --- :mod:`MDAnalysis.coordinates.H5`
===============================================================

This will be a module for reading the H5 format specified by MDTraj:
here is a link:....
as it differs from that of the existing H5MD module.


Classes
-------

.. autoclass:: H5Reader
   :members:
   :inherited-members:

"""

# import json
# import MDAnalysis as mda

import os
import errno
import numpy as np
import struct
import types
import warnings
import logging

from . import base, core
from ..exceptions import NoDataError
from ..due import due, Doi
# from ..lib.util import cached

logger = logging.getLogger("MDAnalysis.coordinates.H5")

try:
    import tables as tb
    # might need simplejson (they use in MDTraj)
except ImportError:
    HAS_TABLES = False

    # Allow building documentation even if h5py is not installed
    # TODO - I don't quite understand how to do this bit here...
#    import types
#
#    class MockH5pyFile:
#        pass
#    h5py = types.ModuleType("h5py")
#    h5py.File = MockH5pyFile

else:
    HAS_TABLES = True


class H5Reader(base.ReaderBase):
    r"""Reader for the H5 - MDTraj format.

    .. versionadded:: 2.2.0

    not implemented:
        delta, skip_timestep as these are not written as far as I know

        see: https://www.mdtraj.org/1.9.8.dev0/api/reporters.html

        mdtraj.reporters.HDF5Reporter(file, reportInterval, coordinates=True, time=True, cell=True, potentialEnergy=True, kineticEnergy=True, temperature=True, velocities=False, atomSubset=None, enforcePeriodicBox=None)

    """

    format = 'H5'
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps'}


    # @due.dcite(MDTraj)
    
    def __init__(self, filename,
                 convert_units=True,
                 dt=None,
                 **kwargs):

        if not HAS_TABLES:
            raise RuntimeError('Please install tables')

        super(H5Reader, self).__init__(filename, convert_units=convert_units, **kwargs)
        self.filename = filename
        self.convert_units = convert_units
        self._cache = dict() # needed for cache decorator

        self.open_trajectory()
        # whcy am I getting an AttributeError trying to set this:
        # self.n_atoms = self.root.coordinates.shape[1]

        # testing stuff: - probably make this into a function called "parse_information" or something...
        self._h5_has_coordinates = False
        self._h5_has_velocites = False
        self._h5_has_time = False
        self._h5_has_cell_angles = False
        self._h5_has_cell_lengths = False
        self._h5_has_topology = False

        # parse file to check for stuff, make sure units are correct, see if there are velocities
        self._parse_file()

        # only want to do dt once:
        # find the time difference between consecutive steps:
        if dt is None:
            dt = np.average(
                self.root.time[1:] - self.root.time[0:-1]
                )
        
        self._ts_kwargs['dt'] = dt

        # initiate the timestep object:
        self.ts = self._Timestep(self.n_atoms, 
                positions=self._h5_has_coordinates,
                velocities=self._h5_has_velocities,
                **self._ts_kwargs)


        self._frame = -1

        self._read_next_timestep(self.ts)
        

        self.ts.dt = dt

        

    def __len__(self):
        '''Get the length of the trajectory'''
        return self.root.coordinates.shape[0]



    def _parse_file(self):

        if not (hasattr(self.root._v_attrs, "application") and self.root._v_attrs.application == "MDTraj"):
            print("Warning, possibly not made by MDTraj") # convert to log message
        else:
            print("yep seems like its from MDTraj")

        if hasattr(self.root, "coordinates"):
            if self.root.coordinates._v_attrs.units != "nanometers":
                raise ValueError("Only works for nanometres at the moment")

            self._h5_has_coordinates = True


        if hasattr(self.root, "velocities"):
            if self.root.velocities._v_attrs.units != "nanometers/picosecond":
                raise ValueError("Only works for nanometres/picosecond at the moment")

            self._h5_has_velocities = True

        if hasattr(self.root, "time"):
            if self.root.time._v_attrs.units != "picoseconds":
                raise ValueError("Only works for picoseconds at the moment")

            self._h5_has_time = True

        if hasattr(self.root, "cell_angles"):
            if self.root.cell_angles._v_attrs.units != "degrees":
                raise ValueError("Only works for degrees at the moment")

            self._h5_has_cell_angles = True

        if hasattr(self.root, "cell_lengths"):
            if self.root.cell_lengths._v_attrs.units != "nanometers":
                raise ValueError("Only works for nanometres at the moment")

            self._h5_has_cell_lengths = True


        # also writes the temperature, potential and kinetic energy - should I check for these?

        if hasattr(self.root, "topology"):
            self._h5_has_topology = True

        
    
    # not sure if this is getting called 
    def close(self):
        """close trajectory"""
        self._file.close()


    def __del__(self):
        # currently this is getting called after the file is closed.
        # for example, when I printed "Calling destructor" here,
        # i get the following...:
            # Closing remaining open files:nvtHiRes.h5...done
            # Calling destructor
        # so maybe the destructor of the tables object is getting called first? 
        # - it appers as though that is the case - is it a problem?
        self.close()


    @staticmethod # might not have to open file if only gets called after __init__?
    def parse_n_atoms(filename, **kwargs):
        with tb.open_file(filename, 'r') as f:
            n_atoms = f.root.coordinates.shape[1] # (nframes, ncoord, dims)
        return n_atoms

    @property
    def n_frames(self):
        """number of frames in a trajectory"""
        return self.root.coordinates.shape[0]

    # should I use @cached like the TXYZ reader?
    @property
    def n_atoms(self):
        """number of atoms in the system"""
        return self.root.coordinates.shape[1]

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)
        """
        return self.ts.dimensions

    @property
    def dt(self):
        """timestep between frames"""
        return self.ts.dt


    def open_trajectory(self):
        # TODO this needs to check the tables version, as the open command is
        # different
        self._file = tb.open_file(self.filename, mode="r")
        self.root = self._file.root


    def _reopen(self):
        """reopen trajectory"""
        self.ts.frame = 0
        self._frame = -1
        self._file.close()
        self.open_trajectory()

    # pass the frame (0 based index) and the ts object:
    def _update_ts(self, frame, ts):
        ts.frame = self._frame # dcd reader has this, not sure if needed?
        ts.time = self.root.time[frame]
        ts.positions = self.root.coordinates[frame]

        if self._h5_has_velocities:
            ts.velocities = self.root.velocities[frame]

        ts.dimensions = np.append(self.root.cell_lengths[frame], self.root.cell_angles[frame])


        if self.convert_units:
            if ts.dimensions is not None:
                self.convert_pos_from_native(ts.dimensions[:3])

            if ts.velocities is not None:
                self.convert_velocities_from_native(ts.velocities)
                
            self.convert_pos_from_native(ts.positions)

            self.convert_time_from_native(ts.time)


        return ts


    def _read_next_timestep(self, ts=None):
        """copy next frame into timestep"""
        if self._frame == self.n_frames - 1:
            raise IOError('trying to go over trajectory limit')
        if ts is None:
            # copy from dcd reader - necessary?
            ts = self.ts.copy()
        self._frame += 1
        ts = self._update_ts(self._frame, ts)
        self.ts = ts
        return ts

    def _read_frame(self, i):
        self._frame = i
        return self._update_ts(i, self.ts)






        
