# opmd2VTK: Converter of openPMD files to VTK

Based on openPMD_viewer and PyVTK

Work in progress: 
- presently only 3D cartesian fields are supported

TBD:
- particles
- CIRC
- units resclaings if needed

### Basic usage :
```python
from opmd_viewer import OpenPMDTimeSeries
from opmd2VTK import Opmd2VTK

ts = OpenPMDTimeSeries('./diags/hdf5/')
conv = Opmd2VTK(ts)

conv.write_fields_vtk(iteration=ts.iterations[-1])
```
