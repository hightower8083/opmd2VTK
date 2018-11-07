# opmd2VTK: Converter of [openPMD](http://www.openpmd.org/#/start) files to  [VTK](https://www.vtk.org) 

The API for OpenPMD-files is based on [openPMD-viewer](https://github.com/openPMD/openPMD-viewer) and possible VTK APIs are [PyVTK](https://github.com/pearu/pyvtk) and 
[tvtk](http://docs.enthought.com/mayavi/tvtk/README.html). PyVTK API is a more lightweight dependance, but tvtk API works with XML format of VTK and is faster.


## Installation
The code was tested only under Linux and MacOS. The recommended installation is through the [Anaconda](https://www.continuum.io/why-anaconda) distribution.
If Anaconda is not your default Python installation, download and install it from [here](https://www.continuum.io/downloads).

It is assumed that openPMD-viewer is already installed. 
- VTK can be installed from anaconda:
```
conda install vtk
```
- for pyvtk API (this step is default and is also done by the installer):
```
pip install pyvtk
```
- tvtk API is a part of Mayavi and can be installed with it:
```
pip install mayavi
```
- clone **opmd2VTK**, `cd` into the folder and install it:
```
git clone https://github.com/hightower8083/opmd2VTK.git
cd opmd2VTK
python setup.py install
```

## Basic usage 

(for pyvtk API replace `opmd2VTK_tvtk` with `opmd2VTK_pyvtk`):

- serial

```python
from opmd_viewer import OpenPMDTimeSeries
from opmd2VTK.opmd2VTK_tvtk import Opmd2VTK

ts = OpenPMDTimeSeries('./diags/hdf5/')
conv = Opmd2VTK(ts)

conv.write_fields_vtk(iteration=ts.iterations[-1])
conv.write_species_vtk(iteration=ts.iterations[-1])
```
- parallel script (needs mpi4py)

```python
from opmd_viewer import OpenPMDTimeSeries
from opmd2VTK.opmd2VTK_tvtk import Opmd2VTK

from mpi4py.MPI import COMM_WORLD  as comm

ts = OpenPMDTimeSeries('./diags/hdf5/')
conv = Opmd2VTK(ts)

for it in ts.iterations[comm.rank::comm.size]:
    conv.write_fields_vtk(iteration=it, zmin_fixed=0)
    conv.write_species_vtk(iteration=it)
```
 
## Useful features :
- Angular resolution of cylindric grid is defined by the argument `Nth` of `write_fields_vtk`;
- `zmin_fixed` allows to fix the origin of visualization domain both for fields and particles (useful for animations);
- argument `CommonMesh=True` groups all fields to the same mesh, while otherwise each component is written to a separate file (resolves the issue of the staggered fields)

## Work in progress

TBD:
- error messages
- debugging

