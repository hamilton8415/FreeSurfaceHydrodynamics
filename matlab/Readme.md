How to compile and use with matlab

In matlab:
>> mex -setup c++
>> mex -I/opt/FreeSurfaceHydrodynamics/include -I/usr/include/eigen3 -L/opt/FreeSurfaceHydrodynamics/lib -lFreeSurfaceHydrodynamics LinearIncidentWave_mex.cpp
