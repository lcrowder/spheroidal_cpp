# Import some kind of testing library
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.0.1)
FetchContent_MakeAvailable(catch)

# Uncomment this to add /extras to the Catch2 path to get access to some cmake macros
#  such as catch_discover_tests
#list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
#include(Catch)

# Need this to ignore a lot of the default tests
set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)
include(CTest)

# Put all testing files here
set( spheroidal_test_list test_simple.cpp
                          test_simple_2.cpp
                          test_simple_3.cpp )

# Add an executable to run all tests
add_executable( all_spheroidal_tests ${spheroidal_test_list} )
target_link_libraries(all_spheroidal_tests spheroidal Catch2::Catch2WithMain fmt::fmt )
add_test( NAME all_spheroidal_tests COMMAND all_spheroidal_tests )
                          
# Add the executable and `add_test` for each individual test file
foreach( test ${spheroidal_test_list} )
    get_filename_component( test_name ${test} NAME_WE )

    add_executable( ${test_name} ${test} )
    target_link_libraries( ${test_name} spheroidal Catch2::Catch2WithMain fmt::fmt )
    
    add_test( NAME ${test_name} COMMAND all_spheroidal_tests ${test_name} )
endforeach()

# If you want to run a specific test from a set, build the file/set it's in, then
#  execute spheroidal_cpp/build/tests/<set name> <test name>.

# TODO: Maybe find a way to print the results of a test to the Doxygen? Unsure how to do that.
# TODO: Figure out the difference between include, find_package, and FetchContent_MakeAvailable