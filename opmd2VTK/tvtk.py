"""
This file is part of opmd2VTK software, which
converts the openPMD files to VTK containers

Copyright (C) 2018, opmd2VTK contributors
Author: Igor Andriyash

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

import numpy as np
import os
from tvtk.api import tvtk as vtk
from tvtk.api import write_data


class Opmd2VTK(opmd2VTKGeneric):
    """
    Main class for the opmd2VTK converter
    Class is initialzed with the OpenPMDTimeSeries object from openPMD-viewer.
    It contains the following public methods:
    - write_fields_vtk
    - write_species_vtk
    and private methods:
    - _convert_mesh_3d
    - _convert_mesh_circ

    For more details, see the corresponding docstrings.
    """

    def write_fields_vtk(self, flds=None, iteration=0,
                         format='binary', zmin_fixed=None,
                         Nth=24, CommonMesh=True):
        """
        Convert the given list of scalar and vector fields from the
        openPMD format to a VTK container, and write it to the disk.

        Parameters
        ----------
        flds: list or None
            List of scalar and vector fields to be converted. If None, it
            converts all available components provided by OpenPMDTimeSeries

        iteration: int
            iteration number to treat (default 0)

        format: str
            format for the VTK file, either 'ascii' or 'binary'

        zmin_fixed: float or None
            When treating the simulation data for the animation, in
            some cases (e.g. with moving window) it is useful to
            fix the origin of the visualization domain. If float number
            is given it will be use as z-origin of the visualization domain

        Nth: int, optional
            Number of nodes of the theta-axis of a cylindric grid in case
            of thetaMode geometry. Note: for high the Nth>>10 the convertion
            may become rather slow.

        CommonMesh: bool
            If True, the same mesh will be used and fields will be converted
            to the scalar and vector types. If False, each component will be
            saved as a separate file with its own grid.
        """
        # Check available fields if comps is not defined
        if flds is None:
            flds = self.ts.avail_fields

        # Register constant parameters
        self.iteration = iteration
        self.zmin_fixed = zmin_fixed
        self.CommonMesh = CommonMesh
        self.Nth = Nth

        # Set grid to None, in order to recomupte it
        self.grid = None

        # Make a numer string for the file to write
        istr = str(self.iteration)
        while len(istr)<7 : istr='0'+istr

        if self.CommonMesh:
            # get and store fields Scalars and Vectors
            vtk_container = []
            for fld in flds:
                field_type = self.ts.fields_metadata[fld]['type']
                if field_type=='vector':
                    fdata, fname = self._convert_field_vec_full(fld)
                    indx = self.grid.point_data.add_array( fdata )
                    self.grid.point_data.get_array(indx).name = fname
                elif field_type=='scalar':
                    fdata, fname = self._convert_field_scl(fld)
                    indx = self.grid.point_data.add_array(fdata)
                    self.grid.point_data.get_array(indx).name = fname

            write_data(self.grid, self.path+'vtk_fields_'+istr)
        else:
            for fld in flds:
                field_type = self.ts.fields_metadata[fld]['type']
                if field_type=='vector':
                    comps = ['x', 'y', 'z']
                    for comp in comps:
                        fld_full = fld + comp
                        file_name = self.path+'vtk_fields_{}_{}'\
                           .format(fld_full, istr)
                        fdata, fname = self._convert_field_vec_comp(fld, comp)
                        indx = self.grid.point_data.add_array(fdata)
                        self.grid.point_data.get_array(indx).name = fname
                        write_data(self.grid, file_name)
                elif field_type=='scalar':
                    file_name = self.path+'vtk_fields_{}_{}'\
                       .format(fld, istr)
                    fdata, fname = self._convert_field_scl(fld)
                    indx = self.grid.point_data.add_array(fdata)
                    self.grid.point_data.get_array(indx).name = fname
                    write_data(self.grid, file_name)

    def write_species_vtk(self, species=None, iteration=0, format='binary',
                          scalars=['ux', 'uy', 'uz', 'w'], select=None,
                          zmin_fixed=None):
        """
        Convert the given list of species from the openPMD format to
        a VTK container, and write it to the disk.

        Parameters
        ----------
        species: list or None
            List of species names to be converted. If None, it
            converts all available species provided by OpenPMDTimeSeries

        scalars: list of strings
            list of values associated with each paricle to be included.
            ex. : 'charge', 'id', 'mass', 'x', 'y', 'z', 'ux', 'uy', 'uz', 'w'

        iteration: int
            iteration number to treat (default 0)

        select: dict
            dictonary to impoer selection on the particles,
            as it is defined in opnePMD_viewer

        format: str
            format for the VTK file, either 'ascii' or 'binary'

        zmin_fixed: float or None
            When treating the simulation data for the animation, in
            some cases (e.g. with moving window) it is useful to
            fix the origin of the visualization domain. If float number
            is given it will be use as z-origin of the visualization domain

        """
        # Check available fields if comps is not defined
        if species is None:
            species = self.ts.avail_species

        # register constants
        self.iteration = iteration
        self.zmin_fixed = zmin_fixed
        self.select = select
        self.scalars = scalars

        # Make a numer string for the file to write
        istr = str(self.iteration)
        while len(istr)<7 : istr='0'+istr
        name_base = self.path+'vtk_specie_{:}_{:}'

        # Convert and save all the species
        for specie in species:
            points, scalars_to_add = self._get_species(specie)

        	# Create the points container
        	self.pts_vtk = vtk.PolyData(points=coords)
	
        	# Create the scalars containers
        	for i, scalar in enumerate(scalars_to_add):
            	indx = self.pts_vtk.point_data.add_array(scalar.astype(self.dtype))
            	self.pts_vtk.point_data.get_array(indx).name = self.scalars[i]
	
            write_data(self.pts_vtk, name_base)

    def _convert_mesh3d(self, origin, resolutions):
        # register the grid VTK container
        self.grid = vtk.ImageData(dimensions=self.dimensions, 
            origin=origin, spacing=resolutions)

    def _convert_mesh_circ(self, dimensions, points):
        # register the grid VTK container
        self.grid = vtk.StructuredGrid(dimensions=dimensions)
        self.grid.points = points
