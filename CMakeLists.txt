cmake_minimum_required(VERSION 3.5)
project(FreeSurfaceHydrodynamics VERSION 1.2.4 LANGUAGES CXX)

configure_file(${PROJECT_SOURCE_DIR}/include/config.h.in
               ${PROJECT_BINARY_DIR}/include/config.h)

#------------------------------------------------------------------------
# Compile as C++17

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(DEFAULT_BUILD_TYPE "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#------------------------------------------------------------------------
# Include cmake
include(${PROJECT_SOURCE_DIR}/cmake/URL.conf.cmake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#--------------------------------------
# Fetch gnuplot-iostream
include(Add_gnuplot-iostream)

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

target_include_directories(${PROJECT_NAME}
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:include>
    PRIVATE
      ${gnuplot-iostream_INCLUDE_DIRS}
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1)


# find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)

# pybind11 method:
pybind11_add_module(fshd src/pybind11/bindings.cpp)
target_link_libraries(fshd PUBLIC FreeSurfaceHydrodynamics)
target_include_directories(fshd
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
)

add_custom_command(TARGET fshd POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:fshd> ${PROJECT_SOURCE_DIR}/fshd_py/fshd
)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)  # In Root, inlcude targets for developers
  add_executable(IncidentWaveExample examples/IncidentWaveExample.cpp)
  target_link_libraries(IncidentWaveExample FreeSurfaceHydrodynamics)
  target_include_directories(IncidentWaveExample PRIVATE ${PROJECT_BINARY_DIR}/include ${gnuplot-iostream_INCLUDE_DIRS})

  add_executable(PlotCoeffsExample examples/PlotCoeffsExample.cpp)
  target_link_libraries(PlotCoeffsExample FreeSurfaceHydrodynamics)
  target_include_directories(PlotCoeffsExample PRIVATE ${PROJECT_BINARY_DIR}/include  ${gnuplot-iostream_INCLUDE_DIRS})

  add_executable(BuoyancyForceExample examples/BuoyancyForceExample.cpp)
  target_link_libraries(BuoyancyForceExample FreeSurfaceHydrodynamics)
  target_include_directories(BuoyancyForceExample PRIVATE ${PROJECT_BINARY_DIR}/include  ${gnuplot-iostream_INCLUDE_DIRS})

  add_executable(GravityForceExample examples/GravityForceExample.cpp)
  target_link_libraries(GravityForceExample FreeSurfaceHydrodynamics)
  target_include_directories(GravityForceExample PRIVATE ${PROJECT_BINARY_DIR}/include  ${gnuplot-iostream_INCLUDE_DIRS})

  add_executable(RadiationForceExample examples/RadiationForceExample.cpp)
  target_link_libraries(RadiationForceExample FreeSurfaceHydrodynamics)
  target_include_directories(RadiationForceExample PRIVATE ${PROJECT_BINARY_DIR}/include  ${gnuplot-iostream_INCLUDE_DIRS})

  add_executable(ExcitingForceExample examples/ExcitingForceExample.cpp)
  target_link_libraries(ExcitingForceExample FreeSurfaceHydrodynamics)
  target_include_directories(ExcitingForceExample PRIVATE ${PROJECT_BINARY_DIR}/include  ${gnuplot-iostream_INCLUDE_DIRS})

  add_executable(MotionExample examples/MotionExample.cpp)
  target_link_libraries(MotionExample FreeSurfaceHydrodynamics)
  target_include_directories(MotionExample PRIVATE ${PROJECT_BINARY_DIR}/include  ${gnuplot-iostream_INCLUDE_DIRS})

  file(COPY examples/example_hydrodynamic_coeffs DESTINATION .)

  install(TARGETS ${PROJECT_NAME} EXPORT FreeSurfaceHydrodynamicsConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})

  install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp")

  install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  install(EXPORT FreeSurfaceHydrodynamicsConfig DESTINATION share/FreeSurfaceHydrodynamics/cmake)

  install(FILES cmake/etc/ld.so.conf.d/${PROJECT_NAME}.conf DESTINATION /etc/ld.so.conf.d/)

  export(TARGETS ${PROJECT_NAME} FILE FreeSurfaceHydrodynamicsConfig.cmake)


  set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

  include(Packing)
endif()
