cmake_minimum_required(VERSION 3.28)

find_package(GTest REQUIRED)

add_executable(MonoChromeWaves_tests MonoChromeWaves.cpp)
target_link_libraries(MonoChromeWaves_tests PRIVATE FreeSurfaceHydrodynamics GTest::gtest_main)

gtest_discover_tests(MonoChromeWaves_tests)
