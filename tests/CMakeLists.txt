# Import some kind of testing library
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.0.1)
FetchContent_MakeAvailable(catch)

# Need to add this to the path to get access to some Catch2 macros
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Need this to ignore a lot of the default tests
set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)
include(CTest)

# Put the filenames without the extension. This is the "set" of tests
set( spheroidal_tests 
    test_simple
    test_simple_2
    test_simple_3
   )

foreach( test ${spheroidal_tests} )
    add_executable( ${test} ${test}.cpp )
    target_link_libraries( ${test} spheroidal Catch2::Catch2WithMain )
    target_include_directories( ${test} PRIVATE ${CMAKE_SOURCE_DIR}/include )
endforeach()

catch_discover_tests( ${spheroidal_tests} ) # Do this instead of a bunch of add_test calls

# If you want to run a specific test from a set, build the file/set it's in, then
#  execute spheroidal_cpp/build/tests/<set name> <test name>.

# TODO: Find a way to properly execute every test all at once without just calling ctest
# TODO: Maybe find a way to print the results of a test to the Doxygen? Unsure how to do that.