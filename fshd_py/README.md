# FreeSurfaceHydrodynamics Python Module (fshd)

## C++ Libs Prerequisite
Please ensure the C++ library has already been built by following the main `README.md`. By default,
the python bindings are built, but they can be turned on/off with the cmake flag
`-DBUILD_PYTHON_BINDINGS=<ON,OFF>`.

Once the C++ library `libFreeSurfaceHydrodynamics.so` and bindings `fshd*.so` are built, the python package,
`fshd_py` can be built.

## Build / Install
Need to grab a few things for python `build --wheel`:
```
$ sudo apt update
$ sudo apt install python3-build python3-wheel python3-venv python3-virtualenv python3-pybind11
```

Navigate to the `fshd_py` python project:
```
$ cd <git_clone_root_dir>/fshd_py
$ python3 -m build --wheel
```

Then, you can either install and use it:
```
$ pip3 install dist/*.whl
$ fshd_motion_example
  -- OR --
import fshd  # in another project
```

Or, you can try it out before installing (this method will also work after installing).
If not installing, you will need to install python prerequisites manually:

```
$ sudo apt update
$ sudo apt install python3-numpy python3-scipy python3-matplotlib
$ pip3 install ode
```

Then (if installed or from the python project folder):
```
$ python3 -m fshd.examples.motion
```

## Examples
All available examples:
```
$ python3 -m fshd.examples.buoyancyforce
$ python3 -m fshd.examples.excitingforce
$ python3 -m fshd.examples.gravityforce
$ python3 -m fshd.examples.incidentwave
$ python3 -m fshd.examples.motion
$ python3 -m fshd.examples.plotcoeffs
$ python3 -m fshd.examples.radiationforce
```

Or, if installed:
```
fshd_buoyancyforce_example
fshd_excitingforce_example
fshd_gravityforce_example
fshd_incidentwave_example
fshd_motion_example
fshd_plotcoeffs_example
fshd_radiationforce_example
```

Check command line options for examples with `-h`:
```
$ fshd_motion_example -h
  -- OR --
$ python3 -m fshd.examples.motion -h

usage: fshd_motion_example [-h] [-a A] [-t T] [-p P] [-b B]

options:
  -h, --help  show this help message and exit
  -a A        (meters) sets the incident wave amplitude
  -t T        (seconds) sets the incident wave period
  -p P        (degrees) sets the incident wave phase angle
  -b B        (degrees) sets the incident wave direction
```
