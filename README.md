This is a library for computing the forces on a body in the presence of waves, and the resulting motions if desired.  The theory is based upon linearied free-surface wave-body intereactions and relies on hydrodynamic coefficients computed by other methods, usually a boundary element method solution performed by a program such as WAMIT.



## How to build
To build the library and examples from source, in the top-level directory of the repository do the following:
   ```
   $ mkdir build
   $ cd build
   $ cmake ..
   $ make
   ```

## Run Examples
There are two example programs that exercise the library in an interactive way.  The source for these is in the examples directory, and after building, they can be run from the build directory as follows:

   ```
   $ IncidentWave_Ex1
   ```
This program exercises the LinearIncidentWave class by computing the wave elevation that represents a mono-chromatic wave specified by randomally selected amplitude, frequency and phase.  The result is plotted at the origin for a range of time, and then also spatially at a selection of time-snapshots.


A more involved examples exercises all of the features of the library, computing the different types of wave-body interaction forces for a floating body, with hydrodynamic behaviors described by included WAMIT output files.  The program has a number of options that can be seen by executing the following in the build directory:
   ```
   $ ./FS_Hydrodynamics_Ex1 -h
   Usage: Test_FS_Hydrodynamics [-cbgremh]
       -c :  Plot Hydrodynamic Coefficients
       -b :  Test Bouyancy Forces
       -g :  Test Gravity Forces
       -r :  Test Radiation Forces
       -e :  Test Exciting Forces
       -m :  Test Motions
       -h :  This help message
   ```

The options perform the following actions:
- '-c':  This reads and plots all of the non-zero hydrodynamic coefficients from the included WAMIT output files for the wave-radiation and wave-exciting forces.  The WAMIT data contains these in both the frequency domain (plotted as a function of frequency), and time-domain (plotted as impulse-response functions).  These hydrodynamic coefficients define the forces computed in the other examples.
- '-b':  This option computes the buoyancy forces in all six degrees of freedom, as a function of displacement that is specified as a sinesoidal motion with time.
- '-g':  This plots the force and moment due to gravity as the body moves in a sinusoidal motion.
- '-r':  This plots all non-zero radiation forces computed from the hydrodynamic coefficients for a prescribed sinesoidal motion in each of six degrees of freedom.  These radiation forces are computed both in the frequency domain (which assumes the body has been oscillating forever), and in the time-domain, which assumes the body has just begun oscillating.  After some small nunber of cycles, the transients at the start of the time-domain simulations should die out and the time-domain and frequency domain results should match.  Slight differences are due to numerical inaccuracies and may be larger or smaller depending on the frequency of the prescribed motion.
- '-e'  This plots six degrees of wave-exciting force for a single frequency and single-direction incoming wave.  Again these are computed in both the time- and frequency-domain, and should match closely.  The transient at the start is not present here because the wave-exciting force computation in the time-domain includes a convolution of the impulse response function over previous times.
- '-m':  This computes the resulting sinesoidal motion resulting from an incoming monochromatic wave by solving the equations of motion for the floating body subject to all of the above forces (buoyancy, gravity, wave-radiation, ane wave-exciting).  Based on the frequency domain hydrodynamic coefficients these results are computed as Response Amplitude Operators that indicate the bodies motion amplitude and phase relative to the amplitude and phase of the incoming wave.  These results are also computed in the time-domain and compared to the frequency domain results. Again these should match after the transients that result from the time-domain initial condition being zero speed at time equal zero. 


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

After installing in either of the above methods, include the following lines in the CMakeLists.txt file:
   ```
   find_package(FreeSurfaceHydrodynamics REQUIRED)
   add_executable(YourProgram YourProgram.cpp)
   target_link_libraries(YourProgram FreeSurfaceHydrodynamics)
   ```


## Theory





## Current limitations


## Documentation
Doxygen generated documentation can be viewed here: https://hamilton8415.github.io/FreeSurfaceHydrodynamics/
