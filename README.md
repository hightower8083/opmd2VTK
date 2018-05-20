# opmd2VTK: Converter of [openPMD](http://www.openpmd.org/#/start) files to  [VTK](https://www.vtk.org) 

Based on [openPMD-viewer](https://github.com/openPMD/openPMD-viewer) and [PyVTK](https://github.com/pearu/pyvtk)

### Work in progress

Done
- 3D cartesian vector and scalar fields
- CIRC cylindrical vector and scalar fields
- particles

TBD:
- error messages
- units rescalings if needed
- other openPMD features, which could have got missed

### Basic usage :
```python
from opmd_viewer import OpenPMDTimeSeries
from opmd2VTK import Opmd2VTK

ts = OpenPMDTimeSeries('./diags/hdf5/')
conv = Opmd2VTK(ts)

conv.write_fields_vtk(iteration=ts.iterations[-1])
conv.write_species_vtk(iteration=ts.iterations[-1])
```
 
### Useful features :
- Angular resolution of cylindric grid is defined by the argument `Nth` of `write_fields_vtk`;
- `zmin_fixed` allows to fix the origin of visualization domain both for fields and particles (useful for animations);
- argument `CommonMesh=True` groups all fields to the same mesh, while otherwise each component is written to a separate file (resolves the issue of the staggered fields)
