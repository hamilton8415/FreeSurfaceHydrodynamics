How to compile and use with matlab

To compile manually and use, in matlab:
>> mex -setup c++
>> mex -I/opt/FreeSurfaceHydrodynamics/include -I/usr/include/eigen3 -L/opt/FreeSurfaceHydrodynamics/lib -lFreeSurfaceHydrodynamics LinearIncidentWave_mex.cpp


Or, this is done automatically by cmake using the find_matlab package.
