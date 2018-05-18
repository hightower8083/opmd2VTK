# opmd2VTK: Converter of openPMD files to VTK

Based on openPMD_viewer and PyVTK

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
