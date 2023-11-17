
1.0.1 - Initial Release 

1.1.0 - February 5, 2023
- Added versioning. CMake Project version now is accessible from classes in the library
- Fixed bug related to not closing windows in some examples.
- Added specification of random seed to Incident Wave constructor
- Added Bretschneider Spectrum.
- Added PM Spectrum options without the unused Tp specification, old methods that had the Tp specification in the prototype remain, but are marked as deprecated and may be removed in the future.
- Added user defined spectrum to LinearIncident Waves.

1.1.1 - February 5, 2023.
- Fixed config.h distrubed in .deb file (was config.h.in, now is config.h as it should be)

1.2.0 - March 1, 2023.
- Changed specification of Incident wave to a shared_ptr, instead of a reference.  API Change.

1.2.1 - March 3, 2023.
- Changed location of installed include files, removed "FreeSurfaceHydrodynamics" prefix directory

1.2.2 - March 5, 2023.
- Fixed missing dependence on pitch/roll for force due to gravity computation

1.3.0 - June 26, 2023.
- Includes python bindings authored by Michael Anderson
- Included randomness in incident wave frequencies to promote longer repeat times
- Changed specification of custom wave-spectrum to be in term of Hz instead of rad/s.

1.3.1 - Sept 32, 2023.
- Changed location of include files to be namespaced as #include "FreeSurfaceHydrodynamics/___.hpp"

