cmake_minimum_required(VERSION 3.14)
project(numerical_propagator_plugin CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to build with Tudat (for full fidelity)
option(WITH_TUDAT "Build with Tudat backend" OFF)

if(WITH_TUDAT)
    add_definitions(-DHAS_TUDAT)
    # TODO: find_package(Tudat) or add_subdirectory for tudat-wasm
endif()

# Plugin library
add_library(numerical_prop STATIC
    src/propagator.cpp
)
target_include_directories(numerical_prop PUBLIC include)

# Tests
add_executable(test_numerical_prop tests/test_propagator.cpp)
target_link_libraries(test_numerical_prop numerical_prop)
target_include_directories(test_numerical_prop PRIVATE include)
