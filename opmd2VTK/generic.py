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

import pyvtk as vtk
import numpy as np
import os

class opmd2VTKGeneric:
    """
    Generic (API-independant) class for the opmd2VTK converter
    Class is initialzed with the OpenPMDTimeSeries object from openPMD-viewer.
    It contains the following private methods:
    - _convert_field_vec_full
    - _convert_field_vec_comp
    - _convert_field_scl
    - _get_opmd_field_3d
    - _get_opmd_field_circ
    - _make_vtk_mesh_3d
    - _get_origin_3d
    - _make_vtk_mesh_circ
    - _get_origin_circ
    - _convert_species

    For more details, see the corresponding docstrings.
    """

    def __init__(self, ts, path_to_dir='./diags/', dtype=np.float32):
        """
        Constuctor of opmd2VTKconverter, based on the OpenPMDTimeSeries
        from openPMD-viewer.

        Parameters
        ----------
        ts: object, OpenPMDTimeSeries
            OpenPMDTimeSeries object initialized with the OpenPMD data

        path_to_dir: string
            The path to the directory where the openPMD files are.

        dtype: numpy dtype
            The data type with which the data should be stored in VTK
        """
        # Register OpenPMDTimeSeries and 3D mesh name
        self.ts = ts
        self.grid = None
        self.zmin_orig = None
        self.dtype = dtype
        self.dimensions = None

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

    def _convert_field_vec_full(self, fld):
        """
        Convert the given scalar or vector fields from the
        openPMD format to a VTK container.

        Parameters
        ----------
        fld: str
            Scalar and vector fields to be converted
            ex. : 'E', 'B', 'J', 'rho'
            converts all available components provided by OpenPMDTimeSeries

        Returns
        -------
        vec: vtk.Vectors
            VTK vector containter
        """
        # Select the tools for the given geometry
        if self.geom=='3dcartesian':
            comps = self.ts.fields_metadata[fld]['axis_labels']
            get_field = self._get_opmd_field_3d
            make_mesh = self._make_vtk_mesh_3d
            get_origin = self._get_origin_3d
        elif self.geom=='thetaMode':
            comps = ['x', 'y', 'z']
            get_field = self._get_opmd_field_circ
            make_mesh = self._make_vtk_mesh_circ
            get_origin = self._get_origin_circ

        flds = []
        origins = []
        for comp in comps:
            fld_data, self.info = get_field(fld, comp=comp)
            flds.append(fld_data)
            make_mesh()
            origins.append(get_origin())

        staggered=False
        for i in range(len(origins)-1):
            if not np.allclose(origins[i+1], origins[i]):
                staggered=True

        if staggered:
            print("Warning: components of {:}-field seem staggered"\
                  .format(fld),
                  "and option CommonMesh=True is used")

        return flds, fld

    def _convert_field_vec_comp(self, fld, comp):
        """
        Convert the given scalar or vector fields from the
        openPMD format to a VTK container.

        Parameters
        ----------
        fld: str
            Scalar and vector fields to be converted
            ex. : 'E', 'B', 'J', 'rho'
            converts all available components provided by OpenPMDTimeSeries

        comp: str
            Field component, for example, 'x', 'y' or 'z'

        Returns
        -------
        scl: vtk.Scalars
            VTK scalar containter
        """
        # Choose the get_field and make_grid finctions for the given geometry
        if self.geom=='3dcartesian':
            get_field = self._get_opmd_field_3d
            make_mesh = self._make_vtk_mesh_3d
        elif self.geom=='thetaMode':
            get_field = self._get_opmd_field_circ
            make_mesh = self._make_vtk_mesh_circ

        fld_data, self.info = get_field(fld, comp=comp)
        make_mesh()

        return fld_data, fld+comp

    def _convert_field_scl(self, fld):
        """
        Convert the given scalar or vector fields from the
        openPMD format to a VTK container.

        Parameters
        ----------
        fld: str
            Scalar and vector fields to be converted
            ex. : 'E', 'B', 'J', 'rho'
            converts all available components provided by OpenPMDTimeSeries

        Returns
        -------
        scl: vtk.Scalars
            VTK scalar containter
        """
        # Choose the get_field and make_grid finctions for the given geometry
        if self.geom=='3dcartesian':
            get_field = self._get_opmd_field_3d
            make_mesh = self._make_vtk_mesh_3d
        elif self.geom=='thetaMode':
            get_field = self._get_opmd_field_circ
            make_mesh = self._make_vtk_mesh_circ

        fld_data, self.info = get_field(fld)
        make_mesh()
        return fld_data, fld

    def _get_opmd_field_3d(self, fld, comp=None):
        """
        Wrapper function to return the 3D array with
        the field.

        Parameters
        ----------
        fld: str
            Scalar and vector fields to be converted
            ex. : 'E', 'B', 'J', 'rho'
            converts all available components provided by OpenPMDTimeSeries

        comp: str or None
            the component of the vector field. If None, assumes the scalar

        Returns
        -------
        fld: numpy.array
            1D NumPy array aligned with the mesh

        info: dict
            fields metadata as returned by the get_field method of
            openPMD TimeSeries
        """
        # Get the particle data
        fld_data, info = self.ts.get_field(fld, coord=comp, slicing=None,
                                           iteration=self.iteration)

        # register the grid dimensions
        self.dimensions = fld_data.shape

        fld_data = fld_data.astype(self.dtype).T.ravel()

        return fld_data, info

    def _get_opmd_field_circ(self, fld, comp=None):
        """
        Wrapper function to return the 3D array with
        the field.

        Parameters
        ----------
        fld: str
            Scalar and vector fields to be converted
            ex. : 'E', 'B', 'J', 'rho'
            converts all available components provided by OpenPMDTimeSeries

        comp: str or None
            the component of the vector field. If None, assumes the scalar

        Returns
        -------
        fld: numpy.array
            1D NumPy array aligned with the mesh

        info: dict
            fields metadata as returned by the get_field method of
            openPMD TimeSeries
        """
        # Load the slice th=0
        fld2d, info = self.ts.get_field(fld, coord=comp,
                                        iteration=self.iteration, theta=0)

        # Get the Z and R axes from the original data
        z, r = info.z, info.r
        r = r[r.size//2:]
        Nz, Nr = z.size, r.size

        # Make Theta axis (half)
        theta = (2*np.pi/self.Nth) * np.arange(self.Nth//2)

        # Allocate the scalar field in the VTK shape (Z,R,Theta)
        fld3d = np.zeros((z.size, r.size, self.Nth+1), dtype=self.dtype)

        # Write the field for th=0
        # Note: array contain a copy of th=0 at th=2*pi, which is
        # needed to close have a closed cylinder
        fld3d[:,:,0] = fld2d[Nr:].T.astype(self.dtype)
        fld3d[:,:,-1] = fld2d[Nr:].T.astype(self.dtype)
        fld3d[:,:,self.Nth//2] = fld2d[:Nr][::-1].T.astype(self.dtype)

        # Write the fields for all other angles
        for i, th in enumerate(theta[1:]):
            fld2d, info = self.ts.get_field(fld, coord=comp, slicing=None,
                                            iteration=self.iteration, theta=th)
            fld3d[:,:,i+1] = fld2d[Nr:].T.astype(self.dtype)
            fld3d[:,:,i+1+self.Nth//2] = fld2d[:Nr][::-1].T.astype(self.dtype)

        fld3d = fld3d.ravel()
        fld2d = None

        return fld3d, info

    def _convert_species(self, species):
        """
        Convert the given species from the openPMD format to a VTK container.

        Parameters
        ----------
        species: str
            Name of the specis
            ex. : 'electrons', 'He+1'

        Returns
        -------
        pts_vtk: vtk.PolyData
            VTK PolyData container with the particles 3D positions

        scalars_vtk: list of vtk.Scalars
            List of scalars associated with the particles
        """
        # Get the particle data
        pts = self.ts.get_particle(var_list=['x', 'y', 'z']+self.scalars,
                                   species=species, iteration=self.iteration,
                                   select=self.select)

        # Split coordinates and scalars
        coords = np.array(pts[:3]).astype(self.dtype).T
        scalars_to_add = pts[3:]

        # Convert microns to meters
        # Note: in further release of openPMD-viewer, coords
        #       are expected to be meters by default
        coords *= 1e-6

        # If zmin_fixed mode is chosen, shift the z-coordinates
        # to match the fields box
        if self.zmin_fixed is not None:
            if self.zmin_orig is None:
                print('zmin_fixed mode can only be use after the fields')
            else:
                coords[:,2] -= self.zmin_orig - self.zmin_fixed

        return coords, scalars_to_add

    def _make_vtk_mesh_3d(self):
        """
        Create a simple 3D mesh using vtk.StructuredPoints method,
        and using the dimensions, origin, resolutions from
        OpenPMDTimeSeries meta-info object (returned as the second
        argument of ts.get_field method).
        """
        # Exit if the grid already exists and CommonMesh is True
        if self.grid is not None and self.CommonMesh:
            return

        # Get origin and resolution of the 3D visualization domain
        origin = self._get_origin_3d()
        resolutions = (self.info.dx, self.info.dy, self.info.dz)
        self._get_mesh_3d(origin, resolutions)

    def _make_vtk_mesh_circ(self):
        """
        Create a cylindric mesh using vtk.StructuredGrid method.
        Note:
            this function operates only once, and all following calls
        """
        # Exit if the grid already exists and CommonMesh is True
        if self.grid is not None and self.CommonMesh:
            return

        # register the z-origin of the grid
        self.zmin_orig = self.info.zmin

        # Get the Z and R axes from the original data
        z, r = self.info.z, self.info.r
        r = r[r.size//2:]
        Nr, Nz = r.size, z.size
        # Make Theta axis (half)
        theta = np.r_[0:2*np.pi:(self.Nth+1)*1j]

        # shift the visualization domain origin if needed
        if self.zmin_fixed is not None:
            z -= z.min() - self.zmin_fixed

        # Make the grid-point (copied from TVTK tutorial)
        points = np.empty([len(theta)*len(r)*len(z),3])
        x_plane = ( np.cos(theta) * r[:,None] ).ravel()
        y_plane = ( np.sin(theta) * r[:,None] ).ravel()
        start = 0
        for iz in range(len(z)):
            end = start + len(x_plane)
            z_plane = z[iz]
            points[start:end,0] = x_plane
            points[start:end,1] = y_plane
            points[start:end,2] = z_plane
            start = end

        # register the grid VTK container
        self._get_mesh_circ(dimensions=(self.Nth+1, Nr, Nz), points)

    def _get_origin_3d(self):
        """
        Get origin of the 3D visualization domain

        Parameters
        ----------
        zmin_fixed: float or None
            When treating the simulation data for the animation, in
            some cases (e.g. with moving window) it is useful to
            fix the origin of the visualization domain. If float number
            is given it will be use as z-origin of the visualization domain

        Returns
        -------
        origin: 3-tuple
            X, Y and Z positions of the visualization  domain
        """
        # register the z-origin of the grid
        self.zmin_orig = self.info.zmin

        # shift the visualization domain origin if needed
        if self.zmin_fixed is None:
            zmin = self.info.zmin
        else:
            zmin = self.zmin_fixed

        xmin = self.info.xmin
        ymin = self.info.ymin

        return (xmin, ymin, zmin)

    def _get_origin_circ(self):
        """
        Get origin of the CIRC visualization domain

        Returns
        -------
        origin: 3-tuple
            X, Y and Z positions of the visualization  domain
        """

        # Get the Z and R axes from the original data
        z, r = self.info.z, self.info.r
        r = r[r.size//2:]

        # shift the visualization domain origin if needed
        if self.zmin_fixed is not None:
            z -= z.min() - self.zmin_fixed

        return (z.min(), -r.max())
