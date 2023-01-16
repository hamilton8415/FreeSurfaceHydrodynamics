cmake_minimum_required(VERSION 3.5)
project(FreeSurfaceHydrodynamics VERSION 1.0.1 LANGUAGES CXX)



set(DEFAULT_BUILD_TYPE "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set(BOOST_USE_STATIC_LIBS OFF)
set(BOOST_USE_MULTITHREADED ON)
set(BOOST_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.71.0 COMPONENTS system iostreams filesystem)

set(CMAKE_INSTALL_PREFIX  /opt/${PROJECT_NAME}/)
include(GNUInstallDirs)

set(LIB_SOURCE_FILES src/FS_Hydrodynamics.cpp src/LinearIncidentWave.cpp src/interp1d.cpp)

add_library(${PROJECT_NAME} SHARED ${LIB_SOURCE_FILES})
target_link_libraries(FreeSurfaceHydrodynamics PUBLIC Eigen3::Eigen Boost::system Boost::iostreams Boost::filesystem)
set_property(TARGET FreeSurfaceHydrodynamics PROPERTY CXX_STANDARD 17)



target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)


set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1)

add_executable(IncidentWaveExample examples/IncidentWaveExample.cpp)
target_link_libraries(IncidentWaveExample FreeSurfaceHydrodynamics)

add_executable(PlotCoeffsExample examples/PlotCoeffsExample.cpp)
target_link_libraries(PlotCoeffsExample FreeSurfaceHydrodynamics)

add_executable(BuoyancyForceExample examples/BuoyancyForceExample.cpp)
target_link_libraries(BuoyancyForceExample FreeSurfaceHydrodynamics)

add_executable(GravityForceExample examples/GravityForceExample.cpp)
target_link_libraries(GravityForceExample FreeSurfaceHydrodynamics)

add_executable(RadiationForceExample examples/RadiationForceExample.cpp)
target_link_libraries(RadiationForceExample FreeSurfaceHydrodynamics)

add_executable(ExcitingForceExample examples/ExcitingForceExample.cpp)
target_link_libraries(ExcitingForceExample FreeSurfaceHydrodynamics)

add_executable(MotionExample examples/MotionExample.cpp)
target_link_libraries(MotionExample FreeSurfaceHydrodynamics)

file(COPY examples/example_hydrodynamic_coeffs DESTINATION .)

install(TARGETS ${PROJECT_NAME} EXPORT FreeSurfaceHydrodynamicsConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

install(EXPORT FreeSurfaceHydrodynamicsConfig DESTINATION share/FreeSurfaceHydrodynamics/cmake)

install(FILES cmake/etc/ld.so.conf.d/${PROJECT_NAME}.conf DESTINATION /etc/ld.so.conf.d/) 

export(TARGETS ${PROJECT_NAME} FILE FreeSurfaceHydrodynamicsConfig.cmake)


set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Packing)

