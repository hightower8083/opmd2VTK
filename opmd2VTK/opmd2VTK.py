"""
This file is part of opmd2VTK software, which
converts the openPMD files to VTK containers

Copyright (C) 2018, opmd2VTK contributors
Authors: Igor Andriyash, ...

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import pyvtk as vtk
import numpy as np
import os

class Opmd2VTK:
    """
    Main class for the opmd2VTK converter
    Class is initialzed with the OpenPMDTimeSeries object from openPMD-viewer.
    It contains the following public methods:
    - write_fields_vtk
    and private methods:
    - _convert_field
    - _get_field_3d
    - _make_grid_3d

    For more details, see the corresponding docstrings.
    """

    def __init__(self, ts, path_to_dir='./diags/'):
        """
        Constuctor of opmd2VTKconverter, based on the OpenPMDTimeSeries
        from openPMD-viewer.

        Parameters
        ----------
        ts: object, OpenPMDTimeSeries
            OpenPMDTimeSeries object initialized with the OpenPMD data

        path_to_dir: string
            The path to the directory where the openPMD files are.
        """
        # Register OpenPMDTimeSeries and 3D mesh name
        self.ts = ts
        self.grid = None

        # Register the geometry type (not supposed to change)
        self.geom = list(self.ts.avail_geom)[0]

        # Check if the /vtk/ folder is in path_to_dir
        # and, if no, created it
        self.path = path_to_dir + 'vtk/'
        if os.path.exists(self.path) == False :
            try:
                os.makedirs(self.path)
            except OSError :
                pass

    def write_fields_vtk(self, comps=None, iteration=0,
                         format='binary', fixed_zmin=None):
        """
        Convert the given list of scalar and vector fields from the
        openPMD format to a VTK container, and write it to the disk.

        Parameters
        ----------
        comps: list or None
            List of scalar and vector fields to be converted. If None, it
            converts all available components provided by OpenPMDTimeSeries

        iteration: int
            iteration number to treat (default 0)

        format: str
            format for the VTK file, either 'ascii' or 'binary'

        fixed_zmin: float or None
            When treating the simulation data for the animation, in
            some cases (e.g. with moving window) it is useful to
            fix the origin of the visualization domain. If float number
            is given it will be use as z-origin of the visualization domain
        """
        # Check available fields if comps is not defined
        if comps is None:
            comps = self.ts.avail_fields

        # Register z-origin of the visualization domain
        self.fixed_zmin = fixed_zmin

        # Convert the fields one by one and store them to the list
        vtk_container = []
        for comp in comps:
            vtk_container.append(self._convert_field(comp,iteration=iteration))

        # Make a numer string for the file to write
        istr = str(iteration)
        while len(istr)<7 : istr='0'+istr

        # Create the VTK data container and write it to the disk
        vtk.VtkData(self.grid, vtk.PointData(*vtk_container))\
            .tofile(self.path+'vtk_fields_'+istr, format=format)

    def _convert_field(self, comp, iteration):
        """
        Convert the given scalar or vector fields from the
        openPMD format to a VTK container.

        Parameters
        ----------
        comp: str
            Scalar and vector fields to be converted
            ex. : 'E', 'B', 'J', 'rho'
            converts all available components provided by OpenPMDTimeSeries

        iteration: int
            iteration number to treat (default 0)
        """
        # Choose the get_field and make_grid finctions for the given geometry
        if self.geom=='3dcartesian':
            get_field = self._get_opmd_field_3d
            make_mesh = self._make_vtk_mesh_3d

        # Converting the vector field
        if self.ts.fields_metadata[comp]['type']=='vector':
            coords = self.ts.fields_metadata[comp]['axis_labels']
            flds = []

            for coord in coords:
                self.fld, self.info = get_field(comp, iteration,
                                                coord=coord)
                flds.append(self.fld.astype(np.float32).T.ravel())
                make_mesh()

            return vtk.Vectors(np.array(flds).T, name=comp)

        # Converting the scalar field
        elif self.ts.fields_metadata[comp]['type']=='scalar':
            self.fld, self.info = get_field(comp, iteration)
            make_mesh()
            return vtk.Scalars(self.fld.astype(np.float32).T.ravel(), name=comp)
        else:
            print('Error')
            return None

    def _get_opmd_field_3d(self, comp, iteration, coord=None):
        """
        Wrapper function to return the 3D array with
        the field.

        Parameters
        ----------
        comp: str
            Scalar and vector fields to be converted
            ex. : 'E', 'B', 'J', 'rho'
            converts all available components provided by OpenPMDTimeSeries

        iteration: int
            iteration number to treat (default 0)

        coord: str or None
            the component of the  vector field. If None, assumes the scalar
        """
        return self.ts.get_field(comp, coord=coord,slicing=None,iteration=iteration)

    def _make_vtk_mesh_3d(self):
        """
        Create a simple 3D mesh using vtk.StructuredPoints method,
        and using the dimensions, origin, resolutions from
        OpenPMDTimeSeries meta-info object (returned as the second
        argument of ts.get_field method).
        Note:
            this function operates only once, and all following calls
        """

        dimensions = self.fld.shape
        if self.fixed_zmin is None:
            origin = (self.info.xmin, self.info.ymin, self.info.zmin)
        else:
            origin = (self.info.xmin, self.info.ymin, fixed_zmin)
        resolutions = (self.info.dx, self.info.dy, self.info.dz)
        self.grid = vtk.StructuredPoints(dimensions, origin, resolutions)
