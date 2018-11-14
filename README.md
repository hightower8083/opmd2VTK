# opmd2VTK: Converter of [openPMD](http://www.openpmd.org/#/start) files to  [VTK](https://www.vtk.org) 

This library provides an easy way to export 3D and quasi-3D simulation data 
from [openPMD](http://www.openpmd.org/#/start) format to [VTK](https://www.vtk.org) 
files. Present API for OpenPMD is based on [openPMD-viewer](https://github.com/openPMD/openPMD-viewer) 
and accepts only HDF5 input format. For VTK, **opmd2VTK** can use either 
[PyVTK](https://github.com/pearu/pyvtk) or 
[tvtk](http://docs.enthought.com/mayavi/tvtk/README.html). 
`PyVTK` is a more lightweight dependancy, while `tvtk` can be much faster as it 
works with a newer VTK format.


## Installation
The code was tested only under Linux and MacOS. The recommended installation is through the [Anaconda](https://www.continuum.io/why-anaconda) distribution.
If Anaconda is not your default Python installation, download and install it from [here](https://www.continuum.io/downloads).


First install the dependencies:
- VTK can be installed from anaconda:
```
conda install vtk
```
- openPMD-viewer can be installed from PyPi:
```
pip install openPMD-viewer
```
- pyvtk can be installed from PyPi:
```
pip install pyvtk
```
- tvtk is a part of Mayavi and can be installed with it from PyPi:
```
pip install mayavi
```


After dependencies are satisfied, clone **opmd2VTK**, `cd` into the folder and install it:
```
git clone https://github.com/hightower8083/opmd2VTK.git
cd opmd2VTK
python setup.py install
```

## Basic usage 

(for pyvtk API replace `opmd2VTK.tvtk` with `opmd2VTK.pyvtk`):

- serial

```python
from opmd_viewer import OpenPMDTimeSeries
from opmd2VTK.tvtk import Opmd2VTK

ts = OpenPMDTimeSeries('./diags/hdf5/')
conv = Opmd2VTK(ts)

conv.write_fields_vtk(iteration=ts.iterations[-1])
conv.write_species_vtk(iteration=ts.iterations[-1])
```
- parallel script (needs mpi4py)

```python
from opmd_viewer import OpenPMDTimeSeries
from opmd2VTK.tvtk import Opmd2VTK

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
- argument `CommonMesh=True` groups all fields to the same mesh, while otherwise each component is written to a separate file (resolves the issue of the staggered fields);
- sampling of field arrays allows in particual reducing the longitudinal resolution, which may often be excess for visualization needs;

## Work in progress

TBD:
- error messages
- debugging

