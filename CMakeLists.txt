### To build from the command line, do the following:
# mkdir build
# cd build
# cmake .. 
# make (<specific target>)
#
### The executable will be in the build directory

cmake_minimum_required(VERSION 3.7...3.26)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(spheroidal_lib VERSION 0.1.0
                       DESCRIPTION "Spheroidal harmonics functions in C++"
                       LANGUAGES C CXX)

# Used to grab fmtlib and Catch2
include(FetchContent) 

# Necessary to link. Consider doing this with FetchContent instead
find_package(GSL REQUIRED)

# Add include directory to the path
include_directories( PRIVATE ${CMAKE_SOURCE_DIR}/include )

# Add source file subdirectories
add_subdirectory(src_gsl_wrapper)
add_subdirectory(src_spheroidal)

# Add example file subdirectories
add_subdirectory(examples)

# Add testing capability.
# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Lets IDEs support folders
    #set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # Docs only available if this is the main app
    find_package(Doxygen)
    if(Doxygen_FOUND)
        add_subdirectory(docs)
    else()
        message(STATUS "Doxygen not found, not building docs")
    endif()

    add_subdirectory(tests)
    enable_testing()
endif()