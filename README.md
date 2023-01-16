This is a library for computing the forces on a body in the presence of waves, and the resulting motions if desired.  The theory is based upon linearied free-surface wave-body intereactions and relies on hydrodynamic coefficients computed by other methods, usually a boundary element method solution performed by a program such as WAMIT.



## How to build
To build the library and examples from source, in the top-level directory of the repository do the following:
   ```
   $ mkdir build
   $ cd build
   $ cmake ..
   $ make
   ```

## Examples
There are a number of example programs that exercise various portions of the library in an interactive way.  The source for these is in the examples directory and can be run from the build directory as indicated below.  Looking at the source code for these examples is a good way to understand the use of the library:


#### Incident Waves
   ```
   $ ./IncidentWaveExamples
   ```
IncidentWaveExamples.cpp exercises the LinearIncidentWave class by computing the wave elevation that represents specified wave.  By default the wave is a monochromatic wave of Amplitude 1.0 meter, Period 12.0 seconds, Phase angle of zero degrees, and direction of 180 degrees (approaching from the positive x direction). These defaults can be over-ridden by the [-a], [-t], [-p], [-b] options respectively.  The '-h' option provides some guidance.  In every case, the resulting wave-elevation is plotted at the origin for a range of time, and then also spatially at a selection of time-snapshots.


#### Hydrodynamic Coefficients 
   ```
   $ ./PlotCoeffsExample
   ```
PlotCoeffsExample.cpp uses the FS_HydroDynamics class to read and plot all of the non-zero hydrodynamic coefficients from the included WAMIT output files for the wave-radiation and wave-exciting forces.  The WAMIT data contains these in both the frequency domain (plotted as a function of frequency), and time-domain (plotted as impulse-response functions).  These hydrodynamic coefficients define the forces computed in the other examples.

#### Buoyancy Forces

   ```
   $ ./BuoyancyForceExample
   ```

BuoyancyForceExample.cpp uses the FS_HydroDynamics class to compute the hydrostatic force and moments  on the buoy in various orientations. the [-a] option can be used to adjust the range of motion.  Options for period and phase do not exist because the resulting and forces and moment are static so the rates of change do not matter. 
Use the [-h] option for more details.


#### Gravity Forces

   ```
   $ ./GravityForceExample
   ```
GravityForceExample.cpp uses the FS_HydroDynamics class to compute the force and moments due to gravity on the buoy in various orientations. the [-a] option can be used to adjust the range of motion.  Options for period and phase do not exist because the resulting and forces and moment are static so the rates of change do not matter.
Use the [-h] option for more details.

#### Wave Radiation Forces

   ```
   $ ./RadiationForceExample
   ```
RadiationForceExample.cpp uses the FS_HydroDynamics class to compute and plot all non-zero radiation forces computed from the hydrodynamic coefficients for a prescribed sinesoidal motion in each of six degrees of freedom.  These radiation forces are computed both in the frequency domain (which assumes the body has been oscillating forever), and also in the time-domain by convolving the motion with the impulse rsponse function.  In the time domain, the body starts from rest and after some small nunber of cycles the transients at the start of the time-domain simulations should die out and the time-domain and frequency domain results should match.  Slight differences are due to numerical inaccuracies and may be larger or smaller depending on the frequency of the prescribed motion.
A default motion amplitude, period, and phase are used, but these can be changed with the [-atp] options.  
Use the [-h] option for more details.

#### Wave Exciting Forces

   ```
   $ ./ExcitingForceExample
   ```
ExcitingForceExample.cpp uses both the FS_HydroDynamics and LinearIncidentWave class to compute and plot wave-exciting force for each of six degrees of freedom for a single frequency and single-direction incoming wave.  Again these are computed in both the time- and frequency-domain, and should match closely.  The transient at the start is not present here because the wave-exciting force computation in the time-domain includes a convolution of the impulse response function over previous times.
A default incident wave amplitude, period, phase, and direction are used, but these can be changed with the [-atpb] options.
Use the [-h] option for more details.



#### Floating Body Motions

   ```
   $ ./MotionExample
   ```
MotionExample.cpp uses both the FS_HydroDynamics and LinearIncidentWave class to compute and plot the resulting sinesoidal motion excited by an incoming monochromatic wave.  This is accomplished by solving the differential equations of motion for the floating body subject to all of the above forces (buoyancy, gravity, wave-radiation, ane wave-exciting).  Based on the frequency domain hydrodynamic coefficients these results are computed as Response Amplitude Operators that indicate the bodies motion amplitude and phase relative to the amplitude and phase of the incoming wave.  These results are also computed in the time-domain and compared to the frequency domain results. Again these should match after the transients that result from the time-domain initial condition being zero speed at time equal zero.
A default incident wave amplitude, period, phase, and direction are used, but these can be changed with the [-atpb] options.
Use the [-h] option for more details.

## How to install
The library can be installed locally after building from source by exectuing the following from the build directory:
  ```
   $ sudo make install
  ``` 

Alternatively, debians of the most recent release are maintained and can be installed after updating apt with the appropriate ppa

  ```
  $ curl -s --compressed "https://hamilton8415.github.io/ppa/KEY.gpg" | gpg --dearmor | sudo tee /etc/apt/trusted.gpg.d/ppa.gpg >/dev/null
  $ sudo curl -s --compressed -o /etc/apt/sources.list.d/my_list_file.list "https://hamilton8415.github.io/ppa/my_list_file.list"
  $ sudo apt update
  $ sudo apt install libfshydrodynamics
  ```

## How to use in another project using CMake

After installing by either of the above methods, include the following lines in the CMakeLists.txt file:
   ```
   find_package(FreeSurfaceHydrodynamics REQUIRED)
   add_executable(YourProgram YourProgram.cpp)
   target_link_libraries(YourProgram FreeSurfaceHydrodynamics)
   ```


## Current limitations
Currently, this library has a number of limitations that may be addressed in future releases.
- Single Wave Direction:  At this time only long-crestline waves from a single-direction can be included.  Adding additional directions through superposition is a straightforward addition, but computational time for wave-exciting forces will grow linearly with the number of directions included.
- Hydrodynamic Coefficient Data:  At this time only WAMIT output files can be read to supply the needed hydrodynamic coefficients.  Further, only the *.1, *.3 files can be read for frequency domain data are supported, and *_IR.1 and *_JR.3 files for the time-domain impulse response function.  Of course, the information between the frequency-domain and time-domain is redundant and it is possible to compute the needed impulse-response functions from the frequency domain data.  This is not done at this time, instead relying on the f2t utility that is part of WAMIT. 
- Axi-symmetric bodies:  This library is currently most useful for axi-symmetric bodies, with this assumption the wave-exciting forces can be computed from a single wave-direction result from the Boundary-Element-Solver (WAMIT), regardless of body orientation in space.  Facilities to read-in hydrodynamic coefficients for more than one wave-direction are not present and would be needed to evaluate the wave-exciting forces on non axi-symmetric bodies of arbitrary orientation.


## Source Code Documentation
Doxygen generated documentation can be viewed here: https://hamilton8415.github.io/FreeSurfaceHydrodynamics/
